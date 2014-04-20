%% aic_henv
% Select the dimension of the envelope subspace using Akaike information
% criterion for the heteroscedastic envelope model.

%% Syntax
%         u = aic_henv(X, Y)
%         u = aic_henv(X, Y, Opts)
%
%% Input
%
% *X*: Group indicators. A matrix with n rows.  X can only have p unique
%  rows, where p is the number of groups. For example, if there 
% are two groups, X can only have 2 different kinds of rows, such as (0, 1)
% and (1, 0), or (1, 0, 10) and (0, 5, 6).  The number of columns is not
% restricted, as long as X only has p unique rows.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% 
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag to print out dimension selection process, 
% Logical 0 or 1. Default value: 0.
% * Opts.table: Flag to tabulate the results, which contains AIC and log
% likelihood for each u. Logical 0 or 1. Default value: 0.
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

r = size(Y, 2);
ic = zeros(r + 1, 1);
llik = zeros(r + 1, 1);

ModelOutput = henv(X, Y, r, Opts);
llik(r + 1) = ModelOutput.l;
ic(r + 1) = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;

for i = 0 : r - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = henv(X, Y, i, Opts);
	llik(i + 1) = ModelOutput.l;
    ic(i + 1) = -2 * ModelOutput.l + 2 * ModelOutput.paramNum;
    
end

[~, u] = min(ic);
u = u - 1;

if tableFlag == 1
    
    fprintf('\n u      log liklihood      AIC\n');
    fprintf('--------------------------------------------\n');
    for i = 0 : r
        fprintf('%2d %15.3f   %12.3f\n', i, llik(i + 1), ic(i + 1));
    end
    fprintf('--------------------------------------------\n');
    
end
