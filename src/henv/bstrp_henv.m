%% bstrp_henv
% Compute bootstrap standard error for the heteroscedastic envelope model. 

%% Syntax
%         bootse = bstrp_henv(X, Y, u, B)
%         bootse = bstrp_henv(X, Y, u, B, Opts)
%
%% Input
%
% *X*: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
% 
% *Y*: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% 
% *u*: Dimension of the envelope subspace.  A positive integer between 0 and
% r.
% 
% *B*: Number of boostrap samples.  A positive integer.
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
%
%% Output
%
% *bootse*: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the heteroscedastic envelope model by bootstrapping the residuals. 

%% Example
%
%         load waterstrider.mat
% 
%         u = lrt_henv(X, Y, 0.01)
%         B = 100;
%         bootse = bstrp_henv(X, Y, u, B)

function bootse = bstrp_henv(X, Y, u, B, Opts)

if nargin < 4
    error('Inputs: X, Y, u and B should be specified!');
elseif nargin == 4
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

DataParameter = make_parameter(X, Y, 'henv');
p = DataParameter.p;
r = DataParameter.r;
n = DataParameter.n;
ng = DataParameter.ng;
ncum = DataParameter.ncum;
ind = DataParameter.ind;

ModelOutput = henv(X, Y, u, Opts);

Yfit = ModelOutput.Yfit;
resi = Y - Yfit;

bootresi = zeros(n, r);
bootBeta = zeros(B, r * p);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of boostrap sample ' int2str(i) '\n']);
    end
    
    for j = 1 : p
        
        if j > 1
            bootresi(ind(ncum(j - 1) + 1 : ncum(j)), :) = ...
            resi(randsample(ind(ncum(j - 1) + 1 : ncum(j)), ng(j), true), :);
        else
            bootresi(ind(1 : ncum(1)), :) = resi(randsample(ind(1 : ng(1)), ng(1), true), :);
        end
        
    end
    
    Yboot = Yfit + bootresi;
    temp = henv(X, Yboot, u, Opts);
    bootBeta(i,:) = reshape(temp.beta, 1, r * p);

end

bootse = reshape(sqrt(diag(cov(bootBeta, 1))), r, p);