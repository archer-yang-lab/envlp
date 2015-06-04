%% aic_xenv
% Use Akaike information criterion to select the dimension of the envelope
% subspace for the reduction on X.

%% Syntax
%         u = aic_xenv(X, Y)
%         u = aic_xenv(X, Y, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations. The predictors must be continuous variables.
% 
% *Y*: Responses. An n by r matrix, r is the number of
% responses. The response can be univariate or multivariate and must be
% continuous variable.
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
% *u*: Dimension of the envelope. An integer between 0 and p.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace for the reduction on X.

%% Example
%
%         load wheatprotein.txt
%         X = wheatprotein(:, 1 : 6);
%         Y = wheatprotein(:, 7);
%         u = aic_xenv(X, Y)

function u = aic_xenv(X, Y, Opts)

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

ModelOutput = xenv(X, Y, p, Opts);
llik(p + 1) = ModelOutput.l;
ic(p + 1) = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;

for i = 0 : p - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = xenv(X, Y, i, Opts);
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
