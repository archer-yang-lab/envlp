%% bstrp_ienv
% Compute bootstrap standard error for the inner envelope model. 

%% Syntax
%         bootse = bstrp_ienv(X, Y, u, B)
%         bootse = bstrp_ienv(X, Y, u, B, Opts)
%
%% Input
%
% *X*: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% 
% *Y*: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% 
% *u*: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.
% 
% *B*: Number of bootstrap samples.  A positive integer.
%
% *Opts*: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out the number of bootstrap samples, 
% logical 0 or 1. Default value: 0.

%% Output
%
% *bootse*: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the inner envelope model by bootstrapping the residuals. 

%% Example
%
%         load irisf.mat
% 
%         u = bic_ienv(X, Y)
%         B = 100;
%         bootse = bstrp_ienv(X, Y, u, B)

function bootse = bstrp_ienv(X, Y, u, B, Opts)

if (nargin < 4)
    error('Inputs: X, Y, B and u should be specified!');
elseif (nargin == 4)
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n r] = size(Y);
p = size(X, 2);

ModelOutput = ienv(X, Y, u, Opts);

Yfit = ones(n, 1) * ModelOutput.alpha' + X * ModelOutput.beta';
resi = Y - Yfit;

bootBeta = zeros(B, r * p);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of boostrap sample ' int2str(i) '\n']);
    end
    
    bootresi = resi(randsample(1 : n, n, true), : );
    Yboot = Yfit + bootresi;
    temp = ienv(X, Yboot, u, Opts);
    bootBeta(i, :) = reshape(temp.beta, 1, r * p);
    
end

bootse = reshape(sqrt(diag(cov(bootBeta, 1))), r, p);