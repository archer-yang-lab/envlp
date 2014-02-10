%% aic_ienv
% Select the dimension of the inner envelope subspace using Akaike
% information criterion.

%% Syntax
%         u = aic_ienv(X, Y)
%         u = aic_ienv(X, Y, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n 
% is the number of observations. The predictors can be univariate or 
% multivariate, discrete or continuous.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.
% 
% *Opts*: A list containing the optional input parameters, to control the
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
% *u*: Dimension of the inner envelope. An integer between 0 and p.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the inner envelope subspace. 

%% Example
%
%         load irisf.mat
%         u = aic_ienv(X, Y)

function u = aic_ienv(X, Y, Opts)

if nargin < 2
    error('Inputs: X, Y should be specified!');
elseif nargin == 2
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

p = size(X, 2);
    
ModelOutput = ienv(X, Y, 0, Opts);
ic = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;
u = 0;


for i = 1 : p
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = ienv(X, Y, i, Opts);
    temp = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;
    
    if temp < ic
        u = i;
        ic = temp;
    end
    
end