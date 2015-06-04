%% sxenv
% Fit the scaled predictor envelope model.

%% Syntax
%         ModelOutput = sxenv(X, Y, u)
%         ModelOutput = sxenv(X, Y, u, Opts)
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
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out number of iterations, logical 0 or 1.
% Default value: 0. 
% * Opts.rep: Number of replicates for scales. This option imposes special 
% structure on scaling parameters. For example, if Opts.rep = [3 4], this 
% means that the first three responses have the same scale and the next 
% four responses share a different scale. The elements of this vector should 
% sum to r. If not specified, the default is [], then all responses will be
% scaled differently. If all responses have the same scale, input [r], then 
% the regular envelope will be applied to the data.
% The input should be a row vector.
% * Opts.Gamma: The initial value for the envelope subspace. A p by u matrix. Default
% value is the one obtained from xenv, with $\Lambda^{-1}X$ being the 
% predictor and Y being the response. 
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
% * ModelOutput.covMatrix: The asymptotic covariance of vec($$\beta$).  A pr by
% pr matrix.  The covariance matrix returned are asymptotic.  For the
% actual standard errors, multiply by 1 / n.
% * ModelOutput.asySE: Asymptotic standard error for elements in $$\beta$ under
% the scaled predictor envelope model.  A p by r matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1 / sqrt(n).
% * ModelOutput.ratio: The asymptotic standard error ratio of the standard
% multivariate linear regression estimator over the scaled predictor envelope
% estimator, for each element in $$\beta$.  A p by r matrix.
% * ModelOutput.paramNum: The number of parameters in the scaled predictor 
% envelope model.  A positive integer.
% * ModelOutput.l: The maximized log likelihood function.  A real number.
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
% # The codes are implemented based on the algorithm in Section 2.2 of Cook 
% and Su (2015).
% # The Grassmann manifold optimization step calls the package sg_min 2.4.3
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).

%% Example
%
% The following codes uses a subset of the chemometrics example in Cook and
% Su (2015) to demonstrate the use of 'sxenv'.
% 
%         load('chemo.mat')
%         X = X(:, [6 11 21 22]);
%         ModelOutput = sxenv(X, Y, 2);
%         ModelOutput.Lambda
%         ModelOutput.ratio
% 
% This example demonstrates the use of Opts.rep.  In this example, the 
% first two predictors are temperatures measured at equal distances along 
% the reactor, and the other two predictors are wall temperature of the
% reactor and the solvent feed rate.
% 
%         load('chemo.mat')
%         X = X(:, [6 11 21 22]);
%         Opts.rep = [2 1 1];
%         ModelOutput = sxenv(X, Y, 2, Opts);
%         ModelOutput.Lambda
%         ModelOutput.ratio

function ModelOutput = sxenv(X, Y, u, Opts)

% Verify and initialize the parameters
%
if nargin < 3
    error('Inputs: X, Y and u should be specified!');
elseif nargin == 3
    Opts = [];
end

X = double(X);
Y = double(Y);

[n, p] = size(X);
[n1, r] = size(Y);

if n ~= n1
    error('The number of observations in X and Y should be equal!');
end

% if p >= r
%     error(['Number of predictors should be less than number of response!' ...
%         '  Please use ordinary least squares.']);
% end

u = floor(u);
if u < 0 || u > p
    error('u should be an integer between [0, p]!');
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

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


%---preparation---
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

if u == 0
        
    eigtem = eig(sigX);
    logDetSigX = sum(log(eigtem(eigtem > 0)));
    eigtem0 = eig(sigY);
    logDetSigY = sum(log(eigtem0(eigtem0 > 0)));
    
    ModelOutput.beta = zeros(p, r);
    ModelOutput.SigmaX = sigX;
    ModelOutput.Lambda = eye(p);
    ModelOutput.Gamma = [];
    ModelOutput.Gamma0 = eye(p);
    ModelOutput.eta = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sigX;
    ModelOutput.muY = mY;
    ModelOutput.muX = mX;
    ModelOutput.SigYcX = sigY;
    ModelOutput.covMatrix = [];
    ModelOutput.asySE = [];
    ModelOutput.ratio = ones(p, r);
    ModelOutput.paramNum = r + p + p * (p + 1) / 2 + r * (r + 1) / 2;  
    ModelOutput.l = - n * (p + r) / 2 * (1 + log(2 * pi)) - n / 2 * (logDetSigX + logDetSigY);
    ModelOutput.n = n;
    
    
elseif u >= (p - (q - 1) / r)
    
    betaOLS = invSigX * sigXY;
    SigYcX = sigY - sigXY' * invSigX * sigXY;
    eigtem = eig(sigX);
    logDetSigX = sum(log(eigtem(eigtem > 0)));
    eigtem0 = eig(SigYcX);
    logDetSigYcX = sum(log(eigtem0(eigtem0 > 0)));
    covMatrix = kron(SigYcX, inv(sigX));
    asyFm = reshape(sqrt(diag(covMatrix)), p, r);
    
    ModelOutput.beta = betaOLS;
    ModelOutput.SigmaX = sigX;
    ModelOutput.Lambda = eye(p);
    ModelOutput.Gamma = eye(p);
    ModelOutput.Gamma0 = [];
    ModelOutput.eta = betaOLS;
    ModelOutput.Omega = sigX;
    ModelOutput.Omega0 = [];
    ModelOutput.muY = mY;
    ModelOutput.muX = mX;
    ModelOutput.SigYcX = SigYcX;
    ModelOutput.covMatrix = covMatrix;
    ModelOutput.asySE = asyFm;
    ModelOutput.ratio = ones(p, r);
    ModelOutput.l = - n * (p + r) / 2 * (1 + log(2 * pi)) - n / 2 * (logDetSigX + logDetSigYcX);
    ModelOutput.paramNum = r + p + r * p + p * (p + 1) / 2 + r * (r + 1) / 2;
    ModelOutput.n = n;
    
else

    maxIter = Opts.maxIter;
	ftol = Opts.ftol;
	gradtol = Opts.gradtol;
	if (Opts.verbose == 0) 
        verbose = 'quiet';
    else
        verbose = 'verbose';
    end
    
    
    if isfield(Opts, 'ite')        
       
        if Opts.ite <= 0
            error('Number of iterations should be a positive integer!');
        end
        
        ite = ceil(Opts.ite);

    else
        
        ite = 10000;
        
    end
    
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
    if isfield(Opts, 'Gamma')
        
        [r2, u2] = size(Opts.Gamma);

        if p ~= r2 || u ~= u2
            error('The size of the initial value of Gamma should be p by u!');
        end
        if rank(Opts.Gamma) < u
            error('The initial value should be full rank!');
        end
        
        Gamma = Opts.Gamma;      

    else
        
        init = xenv(X / Lambda, Y, u, Opts);
        Gamma = init.Gamma;
        
    end
        
    
    for i = 1 : ite
        
        if printFlag == 1
            fprintf(['Current number of iterations ' int2str(i) '\n']);
        end


        DataParameter.Lambda = Lambda;
         
        F = make_F(@F4sxenv, DataParameter);
        dF = make_dF(@dF4sxenv, DataParameter);

        [l, Gamma] = sg_min(F, dF, Gamma, 10, 'prcg', verbose, ftol, gradtol,'quiet');

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
    
    invLambda = diag(1./Ld);
    Gamma0 = grams(nulbasis(Gamma'));
    eta = inv(Gamma' * invLambda * sigX * invLambda * Gamma) * Gamma' * invLambda * sigXY;
    beta = invLambda * Gamma * eta;
    Omega = Gamma' * invLambda * sigX * invLambda * Gamma;
    Omega0 = Gamma0' * invLambda * sigX * invLambda * Gamma0;
    SigmaX = Lambda * (Gamma * Omega * Gamma' + Gamma0 * Omega0 * Gamma0') * Lambda;
    SigYcX = (YC - XC * invLambda * Gamma * eta)' * (YC - XC * invLambda * Gamma * eta) / n;
    eigtem = eig(SigmaX);
    a = sum(log(eigtem(eigtem > 0)));
    b = trace(SigmaX \ sigX);
    eigtem2 = eig(SigYcX);
    c = sum(log(eigtem2(eigtem2 > 0)));
    l = - n * (p + r) * log(2 * pi) / 2 - 0.5 * n * a - 0.5 * n * b - 0.5 * n * c - 0.5 * n * r; 
    
    insigma = inv(SigmaX);
    covMatrix = kron(SigYcX, insigma);
    asyFm = reshape(sqrt(diag(covMatrix)), p, r);
    
    sep1 = p * r;
    J = zeros(p * r + (p + 1) * p / 2, p * r + (p + 1) * p / 2);    
    J(1 : sep1, 1 : sep1) = kron(inv(SigYcX), SigmaX);
    J(sep1 + 1 : end, sep1 + 1 : end) = Expan(p)' * kron(insigma, insigma) * Expan(p) / 2;
    sep2 = p - 1;
    sep3 = p - 1 + r * u;
    sep4 = p - 1 + u * (p + r);
    sep5 = p - 1 + u * (p + r) + u * (u + 1) / 2;
    H = zeros(p * r + (p + 1) * p / 2, p - 1 + r * u + u ^ 2 + p * (p + 1) / 2);
    H(1 : sep1, 1 : sep2) = kron(-eta' * Gamma' * invLambda, inv(Lambda)) * Lmatrix(p);
    H(1 : sep1, sep2 + 1 : sep3) = kron(eye(r), invLambda * Gamma);
    H(1 : sep1, sep3 + 1 : sep4) = kron(eta', inv(Lambda));    
    H(sep1 + 1 : end, 1 : sep2) = 2 * Contr(p) * kron(SigmaX * invLambda, eye(p)) * Lmatrix(p);
    H(sep1 + 1 : end, sep3 + 1 : sep4) = 2 * Contr(p) * (kron(Lambda * Gamma * Omega, Lambda) ...
        - kron(Lambda * Gamma, Lambda * Gamma0 * Omega0 * Gamma0'));
    H(sep1 + 1 : end,  sep4 + 1 : sep5) = Contr(p) * kron(Lambda * Gamma,  Lambda * Gamma) * Expan(u);
    H(sep1 + 1 : end,  sep5 + 1 : end) = Contr(p) * kron(Lambda * Gamma0,  Lambda * Gamma0) * Expan(p - u);
    asyvar = H * pinv(H' * J * H) * H';
    covMatrix = asyvar(1 : p * r, 1 : p * r);
    asySE = reshape(sqrt(diag(covMatrix)), p, r);
    
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
    ModelOutput.covMatrix = covMatrix;
    ModelOutput.asySE = asySE;
    ModelOutput.ratio = asyFm ./ asySE;
    ModelOutput.paramNum = u * r + p * (p + 5) / 2 + r * (r + 3) / 2 - 1;
    ModelOutput.l = l;
    ModelOutput.n = n;
        
end