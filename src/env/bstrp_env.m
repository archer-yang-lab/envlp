%% bstrp_env
% Compute bootstrap standard error for the envelope model. 

%% Syntax
%         bootse = bstrp_env(X, Y, u, B)
%         bootse = bstrp_env(X, Y, u, B, Opts)
%
%% Input
%
% *X*: Predictors, an n by p matrix, p is the number of predictors.  The
% predictors can be univariate or multivariate, discrete or continuous. 
% 
% *Y*: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be
% continuous variables. 
% 
% *u*: Dimension of the envelope subspace.  A positive integer between 0 and
% r.
% 
% *B*: Number of bootstrap samples.  A positive integer.
% 
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out the number of bootstrap samples, 
% logical 0 or 1. Default value: 0.
%
%% Output
%
% *bootse*: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the envelope model by bootstrapping the residuals. 

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1 : 6);
%         alpha = 0.01;
%         u = lrt_env(X, Y, alpha)
%         B = 100;
%         bootse = bstrp_env(X, Y, u, B)

function bootse = bstrp_env(X, Y, u, B, Opts)

if nargin < 4
    error('Inputs: X, Y, B and u should be specified!');
elseif nargin == 4
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n r] = size(Y);
p = size(X, 2);

ModelOutput = env(X, Y, u, Opts);

Yfit = ones(n, 1) * ModelOutput.alpha' + X * ModelOutput.beta';
resi = Y - Yfit;

bootBeta = zeros(B, r * p);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of bootstrap sample ' int2str(i) '\n']);
    end
    
    bootresi = resi(randsample(1 : n, n, true), :);
    Yboot = Yfit + bootresi;
    temp = env(X, Yboot, u, Opts);
    bootBeta(i, :) = reshape(temp.beta, 1, r * p);

end

bootse = reshape(sqrt(diag(cov(bootBeta, 1))), r, p);