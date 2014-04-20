%% aic_penv
% Select the dimension of the partial envelope subspace using Akaike information
% criterion.

%% Syntax
%         u = aic_penv(X, Y)
%         u = aic_penv(X, Y, Opts)
%
%% Input
%
% *X*: A list containing the value of X1 and X2.
% 
% * X.X1: Predictors of main interest. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
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
% * Opts.verbose: Flag to print out dimension selection process. 
% Logical 0 or 1. Default value: 0.
% * Opts.table: Flag to tabulate the results, which contains AIC and log
% likelihood for each u. Logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the partial envelope subspace.  

%% Example
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X.X1 = fiberpaper(:, 7);
%         X.X2 = fiberpaper(:, 5 : 6);
%         u = aic_penv(X, Y)

function u = aic_penv(X, Y, Opts)

if nargin < 2
    error('Inputs: X and Y should be specified!');
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

[~, r] = size(Y);
ic = zeros(r + 1, 1);
llik = zeros(r + 1, 1);

ModelOutput = penv(X, Y, r, Opts);
llik(r + 1) = ModelOutput.l;
ic(r + 1) = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;

for i = 0 : r - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = penv(X, Y, i, Opts);
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
