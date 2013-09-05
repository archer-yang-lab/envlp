%% xenvpls
% Fit the envelope model for the reduction on X using partial least squares
% algorithm.

%% Syntax
%         ModelOutput = xenvpls(X, Y, u)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations. The predictors must be continuous variables.
% 
% *Y*: Responses. An n by r matrix, r is the number of
% responses. The response can be univariate or multivariate and must be
% continuous variable.
%
% *u*: Dimension of the envelope. An integer between 0 and p.
%

%% Output
% 
% *ModelOutput*: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * ModelOutput.beta: The envelope estimator of the regression coefficients $$\beta$. 
% An p by r matrix.
% * ModelOutput.SigX: The envelope estimator of the covariance matrix of X, $$\Sigma_X$.  A p by
% p matrix.
% * ModelOutput.Gamma: The orthogonal basis of the envelope subspace. An p by u
% semi-orthogonal matrix.
% * ModelOutput.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An p by p-u semi-orthogonal matrix.
% * ModelOutput.eta: The coordinates of $$\beta$ with respect to Gamma. An u by r
% matrix.
% * ModelOutput.Omega: The coordinates of $$\Sigma_X$ with respect to Gamma. An u by u
% matrix.
% * ModelOutput.Omega0: The coordinates of $$\Sigma_X$ with respect to Gamma0. An p - u by p - u
% matrix.
% * ModelOutput.mu: The estimated intercept.  An r by 1 vector.
% * ModelOutput.sigYcX: The estimated conditional covariance matrix of Y given X.
% An r by r matrix.
% * ModelOutput.paramNum: The number of parameters in the envelope model.  A positive
% integer.
% * ModelOutput.n: The number of observations in the data.  A positive
% integer.

%% Description
% This function fits the envelope model in the predictor's space,
% by the partial least squares algorithm in Cook et al. (2012). In the
% population level, this algorithm is equivalent to that in xenv.m, which
% uses the maximum likelihood estimation.  In the sample version, the two 
% algorithms are different.  And this algorithm is much faster, which
% provides a root n consistent starting value for the one in xenv.m.

%% Reference
% 
% The codes are implemented based on the algorithm in Section 4.3 of Cook 
% et al (2012).

%% Example
% 
%         load VocabGrowth 
%         X = VocabGrowth(:, 1 : 3); 
%         Y = VocabGrowth(:, 4);
%         m = 5;
%         u = mfoldcv_xenvpls(X, Y, m)
%         ModelOutput = xenvpls(X, Y, u)
%         ModelOutput.beta



function ModelOutput = xenvpls(X, Y, u)

if nargin < 3
    error('Inputs: X, Y and u should be specified!');
end

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

%---preparation---
DataParameter = make_parameter(X, Y, 'xenvpls');

n = DataParameter.n;
p = DataParameter.p;
r = DataParameter.r;
mX = DataParameter.mX;
mY = DataParameter.mY;
sigX = DataParameter.sigX; % Standard estimator of \Sigma_X
sigXY = DataParameter.sigXY;
sigY = DataParameter.sigY;


% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the envelope
% model reduces to the standard model, and it also should be fitted
% differently.


if u > 0 && u < p


    %---Compute \Gamma using partial least squares algorithm---

    Gamma = get_envelope(sigXY, sigX, u);

    %---Compute the rest of the parameters based on \Gamma---
    Gamma0 = grams(nulbasis(Gamma'));
    eta = Gamma' * sigXY;
    Omega = Gamma' * sigX * Gamma;
    Omega0 = Gamma0' * sigX * Gamma0;
    SigX = Gamma * Omega * Gamma' + Gamma0 * Omega0 * Gamma0'; % Envelope estimator of \Sigma_X
    beta = Gamma * inv(Omega) * eta;
    mu = mY - beta' * mX;
    sigYcX = sigY - eta' * inv(Omega) * eta;
    
    ModelOutput.beta = beta;
    ModelOutput.SigX = SigX;
    ModelOutput.Gamma = Gamma;
    ModelOutput.Gamma0 = Gamma0;
    ModelOutput.eta = eta;
    ModelOutput.Omega = Omega;
    ModelOutput.Omega0 = Omega0;
    ModelOutput.mu = mu;
    ModelOutput.sigYcX = sigYcX;
    ModelOutput.paramNum = r + u * r + p * (p + 1) / 2 + r * (r + 1) / 2;
    ModelOutput.n = n;    
    
elseif u == 0
    
    mu = mY; 
    ModelOutput.beta = zeros(p, r);
    ModelOutput.SigX = sigX;
    ModelOutput.Gamma = [];
    ModelOutput.Gamma0 = eye(p);
    ModelOutput.eta = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sigX;
    ModelOutput.mu = mu;
    ModelOutput.sigYcX = sigY;
    ModelOutput.paramNum = r + p * (p + 1) / 2 + r * (r + 1) / 2;
    ModelOutput.n = n;    

elseif u == p
    
    temp = fit_OLS(X, Y);
    beta = temp.betaOLS';
    mu = mY - beta' * mX;
    sigYcX = temp.SigmaOLS;
    
    ModelOutput.beta = beta;
    ModelOutput.SigX = sigX;
    ModelOutput.Gamma = eye(p);
    ModelOutput.Gamma0 = [];
    ModelOutput.eta = beta;
    ModelOutput.Omega = sigX;
    ModelOutput.Omega0 = [];
    ModelOutput.mu = mu;
    ModelOutput.sigYcX = sigYcX;
    ModelOutput.paramNum = r + (p + r) * (p + r + 1) / 2;
    ModelOutput.n = n;    
    
end
    
    
    