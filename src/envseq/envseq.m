%% envseq
% Fit the envelope model using a sequential algorithm.

%% Syntax
%         ModelOutput = envseq(X, Y, u)
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
%% Output
% 
% *ModelOutput*: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * ModelOutput.beta: The envelope estimator of the regression coefficients $$\beta$. 
% An r by p matrix.
% * ModelOutput.Sigma: The envelope estimator of the error covariance matrix.  An r by
% r matrix.
% * ModelOutput.Gamma: The orthogonal basis of the envelope subspace. An r by u
% semi-orthogonal matrix.
% * ModelOutput.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An r by r-u semi-orthogonal matrix.
% * ModelOutput.eta: The coordinates of $$\beta$ with respect to Gamma. An u by p
% matrix.
% * ModelOutput.Omega: The coordinates of Sigma with respect to Gamma. An u by u
% matrix.
% * ModelOutput.Omega0: The coordinates of Sigma with respect to Gamma0. An r-u by r-u
% matrix.
% * ModelOutput.alpha: The estimated intercept in the envelope model.  An r by 1
% vector.
% * ModelOutput.np: The number of parameters in the envelope model.  A positive
% integer.
% * ModelOutput.n: The number of observations in the data.  A positive
% integer.


%% Description
% This function fits the envelope model to the responses and predictors,
% using the maximum likelihood estimation.  When the dimension of the
% envelope is between 1 and r-1, we implemented the algorithm in Cook et
% al. (2010).  When the dimension is r, then the envelope model degenerates
% to the standard multivariate linear regression.  When the dimension is 0,
% it means that X and Y are uncorrelated, and the fitting is different.

%% References
% 
% The codes are implemented based on the sequential algorithm in the lecture notes
% of Cook (2012).

%% Example
% 
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1 : 6);
%         u = mfoldcv_envseq(X, Y, m)
%         ModelOutput = envseq(X, Y, u)
%         ModelOutput.Sigma

function ModelOutput = envseq(X, Y, u)

% Verify the parameters


[n, p] = size(X);
[n1, r] = size(Y);

if n ~= n1
    error('The number of observations in X and Y should be equal!');
end

if p >= r
    error('When the number of responses is less than the number of predictors, the envelope model cannot be applied.');
end

u = floor(u);
if u < 0 || u > r
    error('u should be an integer between [0, r]!');
end


%---preparation---
DataParameter = make_parameter(X, Y, 'env');

n = DataParameter.n;
p = DataParameter.p;
r = DataParameter.r;
mX = DataParameter.mX;
mY = DataParameter.mY;
sigX = DataParameter.sigX;
sigY = DataParameter.sigY;
sigRes = DataParameter.sigRes;
betaOLS = DataParameter.betaOLS;

eigtem = eig(sigY);

F = make_F(@F4env, DataParameter);
dF = make_dF(@dF4env, DataParameter);


% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the envelope
% model reduces to the standard model, and it also should be fitted
% differently.


if u > 0 && u < r

    %---Compute \Gamma using partial least squares algorithm---

    Gamma = get_envelope(betaOLS, sigRes, u);

    %---Compute the rest of the parameters based on \Gamma---
    Gamma0 = grams(nulbasis(Gamma'));
    beta = Gamma * Gamma' * betaOLS;
    alpha = mY - beta * mX;
    eta = Gamma' * beta;
    Omega = Gamma' * sigRes * Gamma;
    Omega0 = Gamma0' * sigY * Gamma0;
    Sigma1 = Gamma * Omega * Gamma';
    Sigma2 = Gamma0 * Omega0 * Gamma0';
    Sigma = Sigma1 + Sigma2;
   
    ModelOutput.beta = beta;
    ModelOutput.Sigma = Sigma;
    ModelOutput.Gamma = Gamma;
    ModelOutput.Gamma0 = Gamma0;
    ModelOutput.eta = eta;
    ModelOutput.Omega = Omega;
    ModelOutput.Omega0 = Omega0;
    ModelOutput.alpha = alpha;
    ModelOutput.np = r + u * p + r * (r + 1) / 2;
    ModelOutput.n = n;
    
elseif u == 0
    
    ModelOutput.beta = zeros(r, p);
    ModelOutput.Sigma = sigY;
    ModelOutput.Gamma = [];
    ModelOutput.Gamma0 = eye(r);
    ModelOutput.eta = [];
    ModelOutput.Omega = [];
    ModelOutput.Omega0 = sigY;
    ModelOutput.alpha = mY;
    ModelOutput.np = r + u * p + r * (r + 1) / 2;
    ModelOutput.n = n;    

elseif u == r
        
    ModelOutput.beta = betaOLS;
    ModelOutput.Sigma = sigRes;
    ModelOutput.Gamma = eye(r);
    ModelOutput.Gamma0 = [];
    ModelOutput.eta = betaOLS;
    ModelOutput.Omega = sigRes;
    ModelOutput.Omega0 = [];
    ModelOutput.alpha = mY - betaOLS * mX;
    ModelOutput.np = r + u * p + r * (r + 1) / 2;
    ModelOutput.n = n;    
    
end
    
    
    