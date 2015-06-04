%% spls
% Implement the scaled SIMPLS algorithm.

%% Syntax
%         ModelOutput = spls(X, Y, u)
%         ModelOutput = spls(X, Y, u, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations.  The predictors must be continuous variables.
%
% *Y*: Responses. An n by r matrix, r is the number of responses. The 
% response can be univariate or multivariate and must be continuous 
% variable.
%
% *u*: Dimension of the envelope. An integer between 0 and p.
%
% *Opts*: A list containing the optional input parameters. If one or several 
% (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out number of iterations, logical 0 or 1.
% Default value: 0. 
% * Opts.Lambda: The initial value for the scales. A p by p diagonal 
% matrix. Default value is the identity matrix. 
%
%% Output
% 
% *ModelOutput*: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * ModelOutput.beta: The scaled predictor envelope estimator of the regression coefficients
% $$\beta$. An p by r matrix.
% * ModelOutput.SigmaX: The scaled envelope estimator of the covariance matrix of X, $$\Sigma_X$.  A p by
% p matrix.
% * ModelOutput.Lambda: The matrix of estimated scales. An p by p diagonal matrix
% with the first diagonal element equal to 1 and other diagonal elements
% being positive.
% * ModelOutput.Gamma: The orthogonal basis of the envelope subspace. A p by u
% semi-orthogonal matrix.
% * ModelOutput.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  A p by p - u semi-orthogonal matrix.
% * ModelOutput.eta: The coordinates of $$\beta$ with respect to Gamma. A u
% by r matrix.
% * ModelOutput.Omega: The coordinates of Sigma with respect to Gamma. A u by u
% matrix.
% * ModelOutput.Omega0: The coordinates of Sigma with respect to Gamma0. A p - u by p - u
% matrix.
% * ModelOutput.muX: The estimated mean of the predictors in the scaled 
% predictor envelope model.  A p by 1 vector.
% * ModelOutput.muY: The estimated mean of the responses in the scaled 
% predictor envelope model.  An r by 1 vector.
% * ModelOutput.sigYcX: The estimated conditional covariance matrix of Y given X.
% An r by r matrix.
% * ModelOutput.paramNum: The number of parameters in the scaled predictor 
% envelope model.  A positive integer.
% * ModelOutput.n: The number of observations in the data.  A positive
% integer.

%% Description
% This function fits the scaled predictor envelope model to the responses 
% and predictors, using the maximum likelihood estimation.  When the dimension of the
% envelope is between 1 and p - 1, we implemented the algorithm in Cook and
% Su (2015).  When the dimension is p, then the scaled predictor envelope model 
% degenerates to the standard multivariate linear regression.  When the
% dimension is 0, it means that X and Y are uncorrelated, and the fitting
% is different.

%% References
% 
% The codes are implemented based on the algorithm in Section 3 of Cook 
% and Su (2015).

%% Example
%
% The following codes uses a subset of the chemometrics example in Cook and
% Su (2015) to demonstrate the use of 'spls'.
% 
%         load('chemo.mat')
%         X = X(:, [6 11 21 22]);
%         ModelOutput = spls(X, Y, 3);
%         ModelOutput.Lambda
% 
% This example demonstrates the use of Opts.rep.  In this example, the 
% first two predictors are temperatures measured at equal distances along 
% the reactor, and the other two predictors are wall temperature of the
% reactor and the solvent feed rate.
% 
%         load('chemo.mat')
%         X = X(:, [6 11 21 22]);
%         Opts.rep = [2 1 1];
%         ModelOutput = spls(X, Y, 3, Opts);
%         ModelOutput.Lambda

function ModelOutput = spls(X, Y, u, Opts)

if nargin < 3
    error('Inputs: X, Y and u should be specified!');
elseif nargin == 3
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

X = double(X);
Y = double(Y);

[n, p] = size(X);
[n1, r] = size(Y);

if n ~= n1
    error('The number of observations in X and Y should be equal!');
end

% if r >= p
%     error('When the number of predictors is less than the number of responses, the envelope model for reduction on X cannot be applied.');
% end

u = floor(u);
if u < 0 || u > p
    error('u should be an integer between [0, p]!');
end

if u > n - 1
    error('The sample size is too small to fit the model with this u.');
end

if isfield(Opts, 'rep')
    ps = sum(Opts.rep);

    if ps ~= p 
        error('The elements in Opts.rep should sum to p');
    end
        
    q = size(Opts.rep, 2);
    if q > p
        error('The numbers of scale parameters cannot exceed the number of predictors');
    end
    rep = Opts.rep;
else
    q = p;
    rep = ones(1, p);
end
cumrep = cumsum(rep, q - 1);


DataParameter = make_parameter(X, Y, 'sxenv');
DataParameter.rep = rep;
n = DataParameter.n;
p = DataParameter.p;
r = DataParameter.r;
XC = DataParameter.XC;
YC = DataParameter.YC;
mX = DataParameter.mX;
mY = DataParameter.mY;
sigX = DataParameter.sigX;
sigY = DataParameter.sigY;
sigXY = DataParameter.sigXY;
sigXcY = DataParameter.sigXcY;
invSigX = DataParameter.invSigX;

if u > 0 && u < p
    
    ite = 10000;
    epsilon = 1e-9; 
    l2 = zeros(1, ite);

    if isfield(Opts, 'Lambda')
        
        [r2, u2] = size(Opts.Lambda);

        if p ~= r2 || p ~= u2
            error('The size of the initial value of Lambda should be p by p!');
        end
        if rank(Opts.Lambda) < p
            error('The initial value of Lambda should be full rank!');
        end
        if sum(sum(abs(Opts.Lambda)))-trace(abs(Opts.Lambda)) ~= 0
            error('The initial value of Lambda should be a diagonal matrix');
        end
        
        d = diag(Opts.Lambda);
        d = d(2 : end)';

    else
        
        d = ones(1, p - 1);
        
    end

    Lambda = diag([1 d]);

    for i = 1 : ite
        
        if printFlag == 1
            fprintf(['Current number of iterations ' int2str(i) '\n']);
        end


        DataParameter.d = d;
        m1 = xenvpls(X / Lambda, Y, u);
        if i == 1
            Gamma = m1.Gamma;
        else
            a = 0.5;
            a = fminsearch(@(a) objfun_spls(a, Gamma, m1.Gamma, DataParameter), a);
            
            Gamma = grams(a * Gamma + (1 - a) * m1.Gamma);
        end
        
        dlam = diag(Lambda);
        d = dlam(cumrep(1 : (q - 1)))';
        [d, l2(i)] = fminsearch(@(d) objfun_sxenv(d, Gamma, DataParameter),  d);

        C = arrayfun(@(x, y) repmat(x, [1 y]), [1 d], rep, 'UniformOutput', false);
        Ld = cell2mat(C);
        Lambda = diag(Ld);
        
        if i > 1 && abs(l2(i) - l2(i - 1)) < epsilon * abs(l2(i))
            break;
        end

    end

    Gamma0 = grams(nulbasis(Gamma'));
    invLambda = diag(1./Ld);
    eta = inv(Gamma' * invLambda * sigX * invLambda * Gamma) * Gamma' * invLambda * sigXY;
    beta = invLambda * Gamma * eta;
    Omega = Gamma' * invLambda * sigX * invLambda * Gamma;
    Omega0 = Gamma0' * invLambda * sigX * invLambda * Gamma0;
    SigmaX = Lambda * (Gamma * Omega * Gamma' + Gamma0 * Omega0 * Gamma0') * Lambda;
    SigYcX = (YC - XC * invLambda * Gamma * eta)' * (YC - XC * invLambda * Gamma * eta) / n;

    ModelOutput.beta = beta;
    ModelOutput.SigmaX = SigmaX;
    ModelOutput.Lambda = Lambda;
    ModelOutput.Gamma = Gamma;
    ModelOutput.Gamma0 = Gamma0;
    ModelOutput.eta = eta;
    ModelOutput.Omega = Omega;
    ModelOutput.Omega0 = Omega0;
    ModelOutput.muY = mY;
    ModelOutput.muX = mX;
    ModelOutput.SigYcX = SigYcX;
    ModelOutput.paramNum = u * r + p * (p + 5) / 2 + r * (r + 3) / 2 - 1;
    ModelOutput.n = n;

elseif u == 0
    
    ModelOutput.beta = zeros(p, r);
    ModelOutput.SigmaX = sigX;
    ModelOutput.Lambda = eye(r);
    ModelOutput.Gamma = [];
    ModelOutput.Gamma0 = eye(r);
    ModelOutput.eta = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sigX;
    ModelOutput.muY = mY;
    ModelOutput.muX = mX;
    ModelOutput.SigYcX = sigY;
    ModelOutput.paramNum = r + p + p * (p + 1) / 2 + r * (r + 1) / 2;  
    ModelOutput.n = n;
    
elseif u == p
    
    betaOLS = invSigX * sigXY;
    SigYcX = sigY - sigXY' * invSigX * sigXY;
    
    ModelOutput.beta = betaOLS;
    ModelOutput.SigmaX = sigX;
    ModelOutput.Lambda = eye(r);
    ModelOutput.Gamma = eye(r);
    ModelOutput.Gamma0 = [];
    ModelOutput.eta = betaOLS;
    ModelOutput.Omega = sigX;
    ModelOutput.Omega0 = [];
    ModelOutput.muY = mY;
    ModelOutput.muX = mX;
    ModelOutput.SigYcX = SigYcX;
    ModelOutput.paramNum = r + p + r * p + p * (p + 1) / 2 + r * (r + 1) / 2;
    ModelOutput.n = n;

end