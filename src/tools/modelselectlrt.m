%% modelselectlrt
% Select the dimension for the envelope family using likelihood ratio
% testing procedure.

%% Syntax
%         u = modelselectlrt(X, Y, alpha, modelType)
%         u = modelselectlrt(X, Y, alpha, modelType, Opts)
% 
%% Input
%
% *X*: Predictors.   The predictors can be univariate or multivariate, 
% discrete or continuous.  
% 
% For model type 'env', 'henv', 'ienv', 
% and 'xenv', X is an n by p matrix, p is the number of
% predictors. 
% 
% For model type 'penv', X is  A list containing the value of X1 and X2.
% 
% * X.X1 (only for 'penv'): Predictors of main interest. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2 (only for 'penv'): Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% *alpha*: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
% 
% *modelType*: A string of characters indicating the model, choices can be 'env',
% 'henv', 'ienv', 'penv' and 'xenv'.
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
% * Opts.table: Flag to tabulate the results, which contains log 
% likelihood, test statistic, degrees of freedom and p-value for each test. 
% Logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the envelope subspace for method 'env', 'henv', 'ienv', 'penv', 
% and 'xenv'.  The likelihood ratio resting procedure does not support
% 'senv', because the scaled envelope models are not nested with the
% standard model.

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1 : 6);
%         alpha = 0.01;
%         modelType = 'env';
%         u = modelselectlrt(X, Y, alpha, modelType)
% 
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X.X1 = fiberpaper(:, 7);
%         X.X2 = fiberpaper(:, 5 : 6);
%         alpha = 0.01;
%         modelType = 'penv';
%         u = modelselectlrt(X, Y, alpha, modelType)

function u = modelselectlrt(X, Y, alpha, modelType, Opts)

% Verify and initialize the parameters
%
if nargin < 4
    error('Inputs: X, Y, alpha and modelType should be specified!');
elseif nargin == 4
    Opts = [];
end

switch(modelType)
    case 'env'
        u = lrt_env(X, Y, alpha, Opts);
    case 'henv'
        u = lrt_henv(X, Y, alpha, Opts);
    case 'ienv'
        u = lrt_ienv(X, Y, alpha, Opts);
    case 'penv'
        u = lrt_penv(X, Y, alpha, Opts);
    case 'xenv'
        u = lrt_xenv(X, Y, alpha, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end