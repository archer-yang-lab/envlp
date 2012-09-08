%% bstrp_envmean
% Compute bootstrap standard error for the envelope estimator of the
% multivariate mean. 

%% Syntax
%         bootse = bstrp_envmean(Y, u, B)
%         bootse = bstrp_envmean(Y, u, B, Opts)
%
%% Input
%
% *Y*: Data matrix. An n by p matrix, p is the dimension of the variable
% and n is number of observations.  
% 
% *u*: Dimension of the envelope subspace.  A positive integer between 0
% and p. 
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
% *bootse*: The standard error for elements in $$\mu$ computed by
% bootstrap.  A p dimensional column vector.

%% Description
% This function computes the bootstrap standard errors for the envelope
% estimator of the multivariate mean by bootstrapping the residuals.  

%% Example
%         load wheatprotein.txt
%         Y = wheatprotein(:, 1 : 6);
%         alpha = 0.01;
%         u = lrt_envmean(Y, alpha)
%         B = 100;
%         bootse = bstrp_envmean(Y, u, B)

function bootse = bstrp_envmean(Y, u, B, Opts)

if nargin < 3
    error('Inputs: Y, B and u should be specified!');
elseif nargin == 3
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

Y = double(Y);

[n p] = size(Y);

ModelOutput = envmean(Y, u, Opts);

Ymean = ones(n, 1) * ModelOutput.mu';
resi = Y - Ymean;

bootBeta = zeros(B, p);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of bootstrap sample ' int2str(i) '\n']);
    end
    
    bootresi = resi(randsample(1 : n, n, true), :);
    Yboot = Ymean + bootresi;
    temp = envmean(Yboot, u, Opts);
    bootBeta(i, :) = temp.mu';

end

bootse = sqrt(diag(cov(bootBeta, 1)));