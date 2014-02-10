%% senv
% Fit the scaled envelope model.

%% Syntax
%         ModelOutput = senv(X, Y, u)
%         ModelOutput = senv(X, Y, u, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p.
%
% *u*: Dimension of the envelope. An integer between 0 and r.
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
%
%% Output
% 
% *ModelOutput*: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * ModelOutput.beta: The scaled envelope estimator of the regression coefficients
% $$\beta$. An r by p matrix.
% * ModelOutput.Sigma: The scaled envelope estimator of the error covariance
% matrix.  An r by r matrix.
% * ModelOutput.Lambda: The matrix of estimated scales. An r by r diagonal matrix
% with the first diagonal element equal to 1 and other diagonal elements
% being positive.
% * ModelOutput.Gamma: The orthogonal basis of the envelope subspace. An r by u
% semi-orthogonal matrix.
% * ModelOutput.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An r by r - u semi-orthogonal matrix.
% * ModelOutput.eta: The coordinates of $$\beta$ with respect to Gamma. A u by p
% matrix.
% * ModelOutput.Omega: The coordinates of Sigma with respect to Gamma. A u by u
% matrix.
% * ModelOutput.Omega0: The coordinates of Sigma with respect to Gamma0. An r - u by r - u
% matrix.
% * ModelOutput.alpha: The estimated intercept in the scaled envelope model.  An r
% by 1 vector.
% * ModelOutput.l: The maximized log likelihood function.  A real number.
% * ModelOutput.covMatrix: The asymptotic covariance of vec($$\beta$).  An rp by
% rp matrix.  The covariance matrix returned are asymptotic.  For the
% actual standard errors, multiply by 1 / n.
% * ModelOutput.asySE: Asymptotic standard error for elements in $$\beta$ under
% the scaled envelope model.  An r by p matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1 / sqrt(n).
% * ModelOutput.ratio: The asymptotic standard error ratio of the standard
% multivariate linear regression estimator over the scaled envelope
% estimator, for each element in $$\beta$.  An r by p matrix.
% * ModelOutput.paramNum: The number of parameters in the scaled envelope model.  A
% positive integer.
% * ModelOutput.n: The number of observations in the data.  A positive
% integer.

%% Description
% This function fits the scaled envelope model to the responses and predictors,
% using the maximum likelihood estimation.  When the dimension of the
% envelope is between 1 and r - 1, we implemented the algorithm in Cook and
% Su (2013).  When the dimension is r, then the scaled envelope model 
% degenerates to the standard multivariate linear regression.  When the
% dimension is 0, it means that X and Y are uncorrelated, and the fitting
% is different.

%% References
% 
% # The codes are implemented based on the algorithm in Section 4.1 of Cook 
% and Su (2013).
% # The Grassmann manifold optimization step calls the package sg_min 2.4.3
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).

%% Example
%
% The following codes produce the results of the test and performance
% example in Cook and Su (2013).
% 
%         load('sales.txt')
%         Y = sales(:, 4 : 7);
%         X = sales(:, 1 : 3);
%         u = bic_env(X, Y)
%         ModelOutput = env(X, Y, u);
%         1 - 1 ./ ModelOutput.ratio
%         u = bic_senv(X, Y)
%         ModelOutput = senv(X, Y, u);
%         ModelOutput.Lambda
%         1 - 1 ./ ModelOutput.ratio
% 
% This example demonstrates the use of opts.rep.  In this example, the 
% first six responses are measured in mg/ml, and the next responses are 
% measured in ug / ml.
% 
%         load('Urine.txt')
%         Y = Urine(:, 3 : 11);   
%         X = Urine(:, 12 : 14); 
%         Opts.rep = [6 3];
%         u = bic_senv(X, Y, Opts)
%         ModelOutput = senv(X, Y, u, Opts)
%         ModelOutput.Lambda


function ModelOutput = senv(X, Y, u, Opts)

% Verify and initialize the parameters
%
if nargin < 3
    error('Inputs: X, Y and u should be specified!');
elseif nargin == 3
    Opts = [];
end

X = double(X);
Y = double(Y);

n = size(X, 1);
[n1, r] = size(Y);

if n ~= n1
    error('The number of observations in X and Y should be equal!');
end

% if p >= r
%     error(['Number of predictors should be less than number of response!' ...
%         '  Please use ordinary least squares.']);
% end

u = floor(u);
if u < 0 || u > r
    error('u should be an integer between [0, r]!');
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

if isfield(Opts, 'rep')
    rs = sum(Opts.rep);

    if rs ~= r 
        error('The elements in Opts.rep should sum to r');
    end
        
    q = size(Opts.rep, 2);
    if q > r
        error('The numbers of scale parameters cannot exceed the number of responses');
    end
    rep = Opts.rep;
else
    q = r;
    rep = ones(1, r);
end

%---preparation---
DataParameter = make_parameter(X, Y, 'senv');
DataParameter.rep = rep;
n = DataParameter.n;
p = DataParameter.p;
r = DataParameter.r;
mX = DataParameter.mX;
mY = DataParameter.mY;
sigX = DataParameter.sigX;
sigY = DataParameter.sigY;
sigRes = DataParameter.sigRes;
betaOLS = DataParameter.betaOLS;
logDetSigY = DataParameter.logDetSigY;

if u == 0
        
    ModelOutput.beta = zeros(r, p);
    ModelOutput.Sigma = sigY;
    ModelOutput.Lambda = eye(r);
    ModelOutput.Gamma = [];
    ModelOutput.Gamma0 = eye(r);
    ModelOutput.eta = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sigY;
    ModelOutput.alpha = mY;
    ModelOutput.l = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * logDetSigY;
    ModelOutput.covMatrix = [];
    ModelOutput.asySE = [];
    ModelOutput.ratio = ones(r, p);
    ModelOutput.paramNum = r + u * p + r * (r + 1) / 2;  
    ModelOutput.n = n;
    
elseif u == r || u >= (q * r - r + 1) / q
    
    covMatrix = kron(inv(sigX), sigRes);
    asyFm = reshape(sqrt(diag(covMatrix)), r, p);
    
    eigtem = eig(sigRes);
    
    ModelOutput.beta = betaOLS;
    ModelOutput.Sigma = sigRes;
    ModelOutput.Lambda = eye(r);
    ModelOutput.Gamma = eye(r);
    ModelOutput.Gamma0 = [];
    ModelOutput.eta = betaOLS;
    ModelOutput.Omega = sigRes;
    ModelOutput.Omega0 = [];
    ModelOutput.alpha = mY - betaOLS * mX;
    ModelOutput.l = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * log(prod(eigtem(eigtem > 0)));
    ModelOutput.covMatrix = covMatrix;
    ModelOutput.asySE = asyFm;
    ModelOutput.ratio = ones(r, p);
    ModelOutput.paramNum = r + r * p + r * (r + 1) / 2;
    ModelOutput.n = n;
    
elseif q == 1
    
    temp = env(X, Y, u);
    
    ModelOutput.beta = temp.beta;
    ModelOutput.Sigma = temp.Sigma;
    ModelOutput.Lambda = eye(r);
    ModelOutput.Gamma = temp.Gamma;
    ModelOutput.Gamma0 = temp.Gamma0;
    ModelOutput.eta = temp.eta;
    ModelOutput.Omega = temp.Omega;
    ModelOutput.Omega0 = temp.Omega0;
    ModelOutput.alpha = temp.alpha;
    ModelOutput.l = temp.l;
    ModelOutput.covMatrix = temp.covMatrix;
    ModelOutput.asySE = temp.asySE;
    ModelOutput.ratio = temp.ratio;
    ModelOutput.paramNum = temp.paramNum;
    ModelOutput.n = temp.n;
    
elseif u > 0 && u < (q * r - r + 1) / q && q > 1

    maxIter = Opts.maxIter;
	ftol = Opts.ftol;
	gradtol = Opts.gradtol;
	if (Opts.verbose == 0) 
        verbose = 'quiet';
    else
        verbose = 'verbose';
    end
    
    ite = 1000;
    epsilon = 1e-9; 
    l2 = zeros(1, ite);
    
    init = env(X, Y, u, Opts);
    Gamma = init.Gamma;
    d = ones(1, q - 1);

    
    
    for i = 1 : ite
        
        if printFlag == 1
            fprintf(['Current number of iterations ' int2str(i) '\n']);
        end
        
        C = arrayfun(@(x, y) repmat(x, [1 y]), [1 d], rep, 'UniformOutput', false);
        Ld = cell2mat(C);

        Lambda = diag(Ld);
        DataParameter.Lambda = Lambda;
        
        F = make_F(@F4senv, DataParameter);
        dF = make_dF(@dF4senv, DataParameter);
        [l, Gamma] = sg_min(F, dF, Gamma, maxIter, 'prcg', verbose, ftol, gradtol);
        
        [d, l2(i)] = fminsearch(@(d) objfun(d, Gamma, DataParameter),  d);
        
        if i > 1 && abs(l2(i) - l2(i - 1)) < epsilon * abs(l2(i))
            break;
        end
        
    end
    
    C = arrayfun(@(x, y) repmat(x, [1 y]), [1 d], rep, 'UniformOutput', false);
    Ld = cell2mat(C);
    Lambda = diag(Ld);
    invLambda = diag(1./Ld);
    Gamma0 = grams(nulbasis(Gamma'));
    eta = Gamma' * invLambda * betaOLS;
    beta = Lambda * Gamma * eta;
    Omega = Gamma' * invLambda * sigRes * invLambda * Gamma;
    Omega0 = Gamma0' * invLambda * sigY * invLambda * Gamma0;
    Sigma = Lambda * (Gamma * Omega * Gamma' + Gamma0 * Omega0 * Gamma0') * Lambda;
        
    ModelOutput.beta = beta;
    ModelOutput.Sigma = Sigma;
    ModelOutput.Lambda = Lambda;
    ModelOutput.Gamma = Gamma;
    ModelOutput.Gamma0 = Gamma0;
    ModelOutput.eta = eta;
    ModelOutput.Omega = Omega;
    ModelOutput.Omega0 = Omega0;
    ModelOutput.alpha = mY - beta * mX;
    ModelOutput.paramNum = init.paramNum + q - 1;
    ModelOutput.l = - 0.5 * l;
    
    %---compute asymptotic variance and get the ratios---
    asyFm = kron(inv(sigX), sigRes);
    asyFm = reshape(sqrt(diag(asyFm)), r, p);
    insigma = inv(Sigma);
    
    sep1 = p * r;
    J = zeros(p * r + (r + 1) * r / 2, p * r + (r + 1) * r / 2);
    J(1 : sep1, 1 : sep1) = kron(sigX, insigma);
    J(sep1 + 1 : end, sep1 + 1 : end) = Expan(r)' * kron(insigma, insigma) * Expan(r) / 2;
    A = zeros(r - 1, q - 1);
    for i = 1 : q - 1
       A(sum(rep(1 : i)), i) = 1; 
    end
    
    sep2 = q - 1;
    sep3 = q - 1 + p * u;
    sep4 = q - 1 + u * (r - u + p);
    sep5 = q - 1 + u * (r - u + p) + u * (u + 1) / 2;
    H = zeros(p * r + (r + 1) * r / 2, q - 1 + p * u + r * (r + 1) / 2);
    
    H(1 : sep1, 1 : sep2) = kron(eta' * Gamma', eye(r)) * Lmatrix(r) * A;
    H(1 : sep1, sep2 + 1 : sep3) = kron(eye(p), Lambda * Gamma);
    H(1 : sep1, sep3 + 1 : sep4) = kron(eta', Lambda * Gamma0);
    
    H(sep1 + 1 : end, 1 : sep2) = Contr(r) * (kron(Lambda * Gamma * Omega * Gamma', eye(r)) ...
        + kron(eye(r), Lambda * Gamma * Omega * Gamma')) * Lmatrix(r) * A ...
        + Contr(r) * (kron(Lambda * Gamma0 * Omega0 * Gamma0', eye(r)) ...
        + kron(eye(r), Lambda * Gamma0 * Omega0 * Gamma0')) * Lmatrix(r) * A;
    H(sep1 + 1 : end, sep3 + 1 : sep4) = 2 * Contr(r) * (kron(Lambda * Gamma * Omega, Lambda * Gamma0) ...
        - kron(Lambda * Gamma, Lambda * Gamma0 * Omega0));
    H(sep1 + 1 : end,  sep4 + 1 : sep5) = Contr(r) * kron(Lambda * Gamma,  Lambda * Gamma) * Expan(u);
    H(sep1 + 1 : end,  sep5 + 1 : end) = Contr(r) * kron(Lambda * Gamma0,  Lambda * Gamma0) * Expan(r - u);
    
    asyvar = H * pinv(H' * J * H) * H';
    covMatrix = asyvar(1 : r * p, 1 : r * p);
    asySE = reshape(sqrt(diag(covMatrix)), r, p);
    
    ModelOutput.covMatrix = covMatrix;
    ModelOutput.asySE = asySE;
    ModelOutput.ratio = asyFm ./ asySE;
    ModelOutput.n = n;
        
end