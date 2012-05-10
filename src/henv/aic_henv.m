%% aic_henv
% Select the dimension of the envelope subspace using Akaike information
% criterion for the heteroscedastic envelope model.

%% Syntax
%         u = aic_henv(X, Y)
%         u = aic_henv(X, Y, Opts)
%
%% Input
%
% *X*: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% 
% *Opts*: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace for the heteroscedastic envelope model. 

%% Example
% 
%         load waterstrider.mat
%         u = aic_henv(X, Y)

function u = aic_henv(X, Y, Opts)

if (nargin < 2)
    error('Inputs: X, Y should be specified!');
elseif (nargin == 2)
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n r] = size(Y);
    
ModelOutput = henv(X, Y, r, Opts);
ic = - 2 * ModelOutput.l + 2 * ModelOutput.np;
u = r;


for i = 0 : r - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = henv(X, Y, i, Opts);
    temp = - 2 * ModelOutput.l + 2 * ModelOutput.np;
    
    if temp < ic
        u = i;
        ic = temp;
    end
    
end
