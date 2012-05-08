%% bstrp_penv
% Compute bootstrap standard error for the partial envelope model. 

%% Syntax
% bootse = bstrp_penv(X, Y, u, B)
% bootse = bstrp_penv(X, Y, u, B, Opts)
%
% Input
%
% X: A list containing the value of X1 and X2.
% 
% * X.X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
% 
% Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% 
% B: Number of boostrap samples.  A positive integer.
% 
% u: Dimension of the partial envelope subspace.  A positive integer between 0 and
% r.
% 
% Opts: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out output, logical 0 or 1. Default value:
% 0.
%
% Output
%
% bootse: The standard error for elements in $$\beta_1$ computed by
% bootstrap.  An r by p1 matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the partial envelope model by bootstrapping the residuals. 

%% Example
% load T7-7.dat
% Y = T7_7(:, 1 : 4);
% Xtemp = T7_7(:, 5 : 7);
% X.X1 = Xtemp(:, 3);
% X.X2 = Xtemp(:, 1 : 2);
% alpha = 0.01;
% u = lrt_penv(X, Y, alpha)
% B = 100;
% bootse = bstrp_penv(X, Y, u, B)

function bootse = bstrp_penv(X, Y, u, B, Opts)

if nargin < 4
    error('Inputs: X, Y, u and B should be specified!');
elseif nargin == 4
    Opts = [];
end

X1 = X.X1;
X2 = X.X2;
[n r] = size(Y);
p1 = size(X1, 2);

ModelOutput = penv(X, Y, u, Opts);

Yfit = ones(n, 1) * ModelOutput.alpha' + X1 * ModelOutput.beta1' + X2 * ModelOutput.beta2';
resi = Y - Yfit;

bootBeta1 = zeros(B, r * p1);

for i = 1 : B
    
    bootresi = resi(randsample(1 : n, n, true), :);
    Yboot = Yfit + bootresi;
    temp = penv(X, Yboot, u, Opts);
    bootBeta1(i, :) = reshape(temp.beta1, 1, r * p1);
    
end

bootse = reshape(sqrt(diag(cov(bootBeta1, 1))), r, p1);