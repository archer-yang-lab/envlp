function ModelOutput = envmean(X, u, Opts)

if nargin < 2
    error('Inputs: X and u should be specified!');
elseif nargin == 2
    Opts = [];
end

[n, p] = size(X);

u = floor(u);
if u < 0 || u > p
    error('u should be an integer between [0, p]!');
end

Opts = make_opts(Opts);

if isfield(Opts, 'init')
    [p2, u2] = size(Opts.init);
    if p ~= p2 || u ~= u2
        error('The size of the initial value should be r by u!');
    end
    if rank(Opts.init) < u2
        error('The initial value should be full rank!');
    end
end


sX = X' * X / n;
mX = mean(X)';
sigX = cov(X, 1);
eigtemX = eig(sX);
logDetSX = log(prod(eigtemX(eigtemX > 0)));

if u > 0 && u < p

    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.sX = sX;
    DataParameter.sigX = sigX;
    DataParameter.logDetSX = logDetSX;

    F = make_F(@F4envmean, DataParameter);
    dF = make_dF(@dF4envmean, DataParameter);

    maxIter = Opts.maxIter;
	ftol = Opts.ftol;
	gradtol = Opts.gradtol;
	if Opts.verbose == 0 
        verbose = 'quiet';
    else
        verbose = 'verbose';
    end
    if ~isfield(Opts, 'init') 
        init = get_Init4envmean(F, u, DataParameter);
    else
        init = Opts.init;
    end
    
    %---Compute \Gamma using sg_min---

    [l Gamma] = sg_min(F, dF, init, maxIter, 'prcg', verbose, ftol, gradtol);

    %---Compute the rest of the parameters based on \Gamma---
    Gamma0 = grams(nulbasis(Gamma'));
    eta = Gamma' * mX;
    mu = Gamma * eta;
    Omega = Gamma' * sigX * Gamma;
    Omega0 = Gamma0' * sX * Gamma0;
    Sigma = Gamma * Omega * Gamma' + Gamma0 * Omega0 * Gamma0';

    %---compute asymptotic variance and get the ratios---
    asyFm = sqrt(diag(Sigma));
    temp = kron(eta * eta', inv(Omega0)) + kron(Omega, inv(Omega0))... 
        + kron(inv(Omega), Omega0) - 2 * kron(eye(u), eye(p - u));  
    covMatrix = Gamma * Omega * Gamma' + kron(eta', Gamma0) * inv(temp) * kron(eta, Gamma0');
    asyEnv = sqrt(diag(covMatrix));

    ModelOutput.mu = mu;
    ModelOutput.Sigma = Sigma;
    ModelOutput.Gamma = Gamma;
    ModelOutput.Gamma0 = Gamma0;
    ModelOutput.eta = eta;
    ModelOutput.Omega = Omega;
    ModelOutput.Omega0 = Omega0;
    ModelOutput.l = - 0.5 * l;
    ModelOutput.covMatrix = covMatrix;
    ModelOutput.asyEnv = asyEnv;
    ModelOutput.ratio = asyFm;
    ModelOutput.np = u + p * (p + 1) / 2;
    ModelOutput.n = n;
    
elseif u == 0
    
    ModelOutput.mu = zeros(p, 1);
    ModelOutput.Sigma = sX;
    ModelOutput.Gamma = [];
    ModelOutput.Gamma0 = eye(p);
    ModelOutput.eta = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sX;
    ModelOutput.l = - n * p / 2 * (1 + log(2 * pi)) - n / 2 * logDetSX;
    ModelOutput.covMatrix = [];
    ModelOutput.asyEnv = [];
    ModelOutput.ratio = ones(p, 1);
    ModelOutput.np = u + p * (p + 1) / 2;
    ModelOutput.n = n;  
    
elseif u == p
    
    eigtem = eig(sigX);
    logDetSigX = log(prod(eigtem(eigtem > 0)));
    
    ModelOutput.mu = mX;
    ModelOutput.Sigma = sigX;
    ModelOutput.Gamma = eye(p);
    ModelOutput.Gamma0 = [];
    ModelOutput.eta = mX;
    ModelOutput.Omega = sigX;
    ModelOutput.Omega0 = [];
    ModelOutput.l = - n * p / 2 * (1 + log(2 * pi)) - n / 2 * logDetSigX;
    ModelOutput.covMatrix = sigX;
    ModelOutput.asyEnv = sqrt(diag(sigX));
    ModelOutput.ratio = ones(p, 1);
    ModelOutput.np = u + p * (p + 1) / 2;
    ModelOutput.n = n; 
    
end