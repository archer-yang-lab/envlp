%% ienv
% Fit the inner envelope model.

%% Syntax
%         ModelOutput = ienv(X, Y, u)
%         ModelOutput = ienv(X, Y, u, Opts)
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
% *u*: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.
%
% *Opts*: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out Grassmann manifold optimization 
% process, logical 0 or 1. Default value: 0.
%
%% Output
% 
% *ModelOutput*: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * ModelOutput.beta: The envelope estimator of the regression coefficients $$\beta$. 
% An r by p matrix.
% * ModelOutput.Sigma: The envelope estimator of the error covariance matrix.  An r by
% r matrix.
% * ModelOutput.Gamma1: The orthogonal basis of the inner envelope subspace. An r by u
% semi-orthogonal matrix.
% * ModelOutput.Gamma0: The orthogonal basis of the complement of the inner envelope
% subspace.  An r by r-u semi-orthogonal matrix.
% * ModelOutput.eta1: The transpose of the coordinates of $$\beta$ with respect to
% Gamma1. An p by u matrix.
% * ModelOutput.B: An (r - u) by (p - u) semi-orthogonal matrix, so that (Gamma,
% Gamma0 * B) spans $$\beta$.
% * ModelOutput.eta2: The transpose of the coordinates of $$\beta$ with respect to
% Gamma0. An p by (p-u) matrix.
% * ModelOutput.Omega1: The coordinates of Sigma with respect to Gamma1. An u by u
% matrix.
% * ModelOutput.Omega0: The coordinates of Sigma with respect to Gamma0. An r - u by r - u
% matrix.
% * ModelOutput.alpha: The estimated intercept in the inner envelope model.  An r by 1
% vector.
% * ModelOutput.l: The maximized log likelihood function.  A real number.
% * ModelOutput.covMatrix: The asymptotic covariance of vec($$\beta$).  An rp by
% rp matrix.  The covariance matrix returned are asymptotic.  For the
% actual standard errors, multiply by 1/n.
% * ModelOutput.asyIenv: Asymptotic standard error for elements in $$\beta$ under
% the inner envelope model.  An r by p matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1 / sqrt(n).
% * ModelOutput.ratio: The asymptotic standard error ratio of the standard multivariate 
% linear regression estimator over the inner envelope estimator, for each element 
% in $$\beta$.  An r by p matrix.
% * ModelOutput.np: The number of parameters in the inner envelope model.  A positive
% integer.
% * ModelOutput.n: The number of observations in the data.  A positive
% integer.

%% Description
% This function fits the inner envelope model to the responses and predictors,
% using the maximum likelihood estimation.  When the dimension of the
% envelope is between 1 and p-1, we implemented the algorithm in Su and
% Cook (2012).  When the dimension is p, then the inner envelope model degenerates
% to the standard multivariate linear regression.  When the dimension is 0,
% it means that X and Y are uncorrelated, and the fitting is different.

%% References
% 
% # The codes is implemented based on the algorithm in Su and Cook (2012).
% # The Grassmann manifold optimization step calls the package sg_min 2.4.1
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).

%% Example
%
% The following codes gives the results of the Fisher's iris data example
% in Su and Cook (2012).
% 
%         load irisf.mat
% 
%         u = bic_env(X, Y)
%         d = bic_ienv(X, Y)
%         ModelOutput = ienv(X, Y, d)
%         1 - 1 ./ ModelOutput.ratio

function ModelOutput = ienv(X, Y, u, Opts)

if (nargin < 3)
    error('Inputs: X, Y and u should be specified!');
elseif (nargin == 3)
    Opts = [];
end

[n,p] = size(X);
[n1,r] = size(Y);

if (n ~= n1)
    error('The number of observations in X and Y should be equal!');
end

if (p >= r)
    error('When the number of responses is less than the number of predictors, the inner envelope model cannot be applied.');
end

u = floor(u);
if (u < 0 || u > p)
    error('u should be an integer in [0, p]!');
end

Opts = make_opts(Opts);

if isfield(Opts, 'init')
    [r2,u2] = size(Opts.init);

    if (r ~= r2 || u ~= u2)
        error('The size of the initial value should be r by u!');
    end

    if (rank(Opts.init) < u2)
        error('The initial value should be full rank!');
    end
end

DataParameter = make_parameter(X, Y, 'ienv');
n = DataParameter.n;
p = DataParameter.p;
r = DataParameter.r;
XC = DataParameter.XC;
YC = DataParameter.YC;
mX = DataParameter.mX;
mY = DataParameter.mY;
sigX = DataParameter.sigX;
sigY = DataParameter.sigY;
sigFit = DataParameter.sigFit;
sigRes = DataParameter.sigRes;
betaOLS = DataParameter.betaOLS;

if u == p
    
    temp = env(X, Y, p, Opts);
    ModelOutput.beta = temp.beta;
    ModelOutput.Sigma = temp.Sigma;
    ModelOutput.Gamma1 = temp.Gamma;
    ModelOutput.Gamma0 = temp.Gamma0;
    ModelOutput.B = [];
    ModelOutput.eta1 = temp.eta';
    ModelOutput.eta2 = [];
    ModelOutput.Omega1 = temp.Omega;
    ModelOutput.Omega0 = temp.Omega0;
    ModelOutput.alpha = temp.alpha;
    ModelOutput.l = temp.l;
    ModelOutput.np = temp.np;
    ModelOutput.covMatrix = temp.covMatrix;
    ModelOutput.asyIenv = temp.asyEnv;
    ModelOutput.ratio = temp.ratio;
    ModelOutput.n = n;
    
elseif u == 0
    
    eigtem = eig(sigY);
    ModelOutput.beta = zeros(r, p);
    ModelOutput.Sigma = sigY;
    ModelOutput.Gamma1 = [];
    ModelOutput.Gamma0 = eye(r);
    ModelOutput.B = eye(r);
    ModelOutput.eta1 = [];
    ModelOutput.eta2 = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sigY;
    ModelOutput.alpha = mY;
    ModelOutput.l = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * log(prod(eigtem(eigtem > 0)));
    ModelOutput.covMatrix = [];
    ModelOutput.asyIenv = [];
    ModelOutput.ratio = ones(r, p);
    ModelOutput.np = r + u * p + r * (r + 1) / 2;  
    ModelOutput.n = n;

else 
    
    eigtem = eig(sigRes);
    
    F = make_F(@F4ienv, DataParameter);
    dF = make_dF(@dF4ienv, DataParameter);

    maxIter = Opts.maxIter;
	ftol = Opts.ftol;
	gradtol = Opts.gradtol;
	if (Opts.verbose == 0) 
        verbose = 'quiet';
    else
        verbose = 'verbose';
    end
    if ~isfield(Opts, 'init') 
        init = get_Init(F, X, Y, u, DataParameter);
    else
        init = Opts.init;
    end
    
    [l Gamma1] = sg_min(F, dF, init, maxIter, 'prcg', verbose, ftol, gradtol);
    
    Gamma0 = grams(nulbasis(Gamma1'));
    eta1 = (Gamma1' * betaOLS)';
    Omega1 = (YC * Gamma1 - XC * eta1)' * (YC * Gamma1 - XC * eta1) / n;
    
    [Vtmp Dtmp] = eig(inv(Gamma0' * sigRes * Gamma0));
    temp1 = Vtmp * diag(sqrt(diag(Dtmp))) * Vtmp';
    temp2 = Gamma0' * sigFit * Gamma0;
    temp3 = Vtmp * diag(diag(Dtmp) .^ ( - 0.5)) * Vtmp';
    [Vt Lambdat] = eig(temp1 * temp2 * temp1);
    [Lambda ind] = sort(diag(Lambdat), 'descend');
    V = Vt(:, ind);
    K = diag([zeros(1, p - u), Lambda(p - u + 1 : end)']);
    Omega0 = Gamma0' * sigRes * Gamma0 + temp3 * V * K * V' * temp3;
    
    [Vtmp Dtmp] = eig(inv(Omega0));
    Omega0inv12 = Vtmp * diag(sqrt(diag(Dtmp))) * Vtmp;
    [Vtmp Dtmp] = eig(Omega0inv12 * Gamma0' * sigFit * Gamma0 * Omega0inv12);
    [Ds ind] = sort(Dtmp, 'descend');
    B = grams(inv(Omega0inv12) * Vtmp(:, ind(1 : p - u)));
    eta2 = (inv(B' * inv(Omega0) * B) * B' * inv(Omega0) * Gamma0' * betaOLS)';
    beta = Gamma1 * eta1' + Gamma0 * B * eta2';
    Sigma = Gamma1 * Omega1 * Gamma1' + Gamma0 * Omega0 * Gamma0';
    alpha = mY - beta * mX;
    
    ModelOutput.beta = beta;
    ModelOutput.Sigma = Sigma;
    ModelOutput.Gamma1 = Gamma1;
    ModelOutput.Gamma0 = Gamma0;
    ModelOutput.B = B;
    ModelOutput.eta1 = eta1;
    ModelOutput.eta2 = eta2;
    ModelOutput.Omega1 = Omega1;
    ModelOutput.Omega0 = Omega0;
    ModelOutput.alpha = alpha;
    ModelOutput.np = p ^ 2 + (p - u) * (r - p) + r * (r + 1) / 2;
    ModelOutput.l = - n * r / 2 * (1 + log(2 * pi)) - n / 2 * log(prod(eigtem(eigtem > 0))) - n / 2 * l;
    
    
    %-----Compute asymptotic variance for inner envelope model-----
    %-----Standard Model-----
    asyFm = kron(inv(sigX), Sigma);
    asyFm = reshape(sqrt(diag(asyFm)), r, p);
    
    %-----Inner Envelope Model-----
    B0 = grams(nulbasis(B'));

    H = zeros(p * r + (r + 1) * r / 2, p * p + (p - u) * (r - p) + u * r + u * (u + 1) / 2 + (r - u) * (r - u + 1) / 2);

    H(1 : p * r, 1 : p * u) = kron(eye(p), Gamma1);
    H(1 : p * r, p * u + 1 : p * p) = kron(eye(p), Gamma0 * B);
    H(1 : p * r, p * p + 1 : p * p + (p - u) * (r - p)) = kron(eta2, Gamma0 * B0);
    H(1 : p * r, p * p + (p - u) * (r - p) + 1 : p * p + (p - u) * (r - p) + r * u) ...
		= kron(eta1, eye(r)) - kron(eta2 * B' * Gamma0', Gamma1) * Kpd(r, u);
    H(p * r + 1 : end, p * p + (p - u) * (r - p) + 1 : p * p + (p - u) * (r - p) + r * u) ...
		= 2 * Contr(r) * (kron(Gamma1 * Omega1, eye(r)) - kron(Gamma1, Gamma0 * Omega0 * Gamma0'));
    H(p * r + 1 : end, p * p + (p - u) * (r - p) + r * u + 1 : p * p + (p - u) * (r - p) + r * u + u * (u + 1) / 2) ...
		= Contr(r) * kron(Gamma1, Gamma1) * Expan(u);
    H(p * r + 1 : end, p * p + (p - u) * (r - p) + r * u + u * (u + 1) / 2 + 1 : end) ... 
		= Contr(r) * kron(Gamma0, Gamma0) * Expan(r - u);

    insigma = inv(Sigma);
    J = zeros(p * r + (r + 1) * r / 2, p * r + (r + 1) * r / 2);
    J(1 : p * r, 1 : p * r) = kron(sigX, insigma);
    J(p * r + 1 : end, p * r + 1 : end) = Expan(r)' * kron(insigma, insigma) * Expan(r) / 2;
    asyv = H * pinv(H' * J * H) * H';
    covMatrix = asyv(1 : r * p, 1 : r * p);
    asyIenv = reshape(sqrt(diag(covMatrix)), r, p);    
    
    ModelOutput.covMatrix = covMatrix;
    ModelOutput.asyIenv = asyIenv;
    ModelOutput.ratio = asyFm ./ asyIenv;
    ModelOutput.n = n;
    
end


