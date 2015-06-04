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
% * Opts.verbose: Flag to print out dimension selection process. 
% Logical 0 or 1. Default value: 0.
% * Opts.table: Flag to tabulate the results, which contains AIC and log
% likelihood for each u. Logical 0 or 1. Default value: 0.
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

if isfield(Opts, 'table')
    if (Opts.table ~= 1)
        tableFlag = 0;
    else
        tableFlag = 1;
    end
else
    tableFlag = 0;
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

p = size(X, 2);
ic = zeros(p + 1, 1);
llik = zeros(p + 1, 1);

ModelOutput = ienv(X, Y, 0, Opts);
llik(1) = ModelOutput.l;
ic(1) = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;

for i = 1 : p
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = ienv(X, Y, i, Opts);
	llik(i + 1) = ModelOutput.l;
    ic(i + 1) = -2 * ModelOutput.l + 2 * ModelOutput.paramNum;
    
end

[~, u] = min(ic);
u = u - 1;

if tableFlag == 1
    
    fprintf('\n u      log likelihood      AIC\n');
    fprintf('--------------------------------------------\n');
    for i = 0 : p
        fprintf('%2d %15.3f   %12.3f\n', i, llik(i + 1), ic(i + 1));
    end
    fprintf('--------------------------------------------\n');
    
end