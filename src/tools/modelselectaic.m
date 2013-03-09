%% modelselectaic
% Select the dimension for the envelope family using Akaike information
% criteria.

%% Syntax
%         u = modelselectaic(X, Y, modelType)
%         u = modelselectaic(X, Y, modelType, Opts)
% 
%% Input
%
% *X*: Predictors.   The predictors can be univariate or multivariate, 
% discrete or continuous.  
% 
% For model type for method 'env', 'henv', 'ienv',
% ' senv', and 'xenv'. X is an n by p matrix, p is the number of
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
% *modelType*: A string of characters indicting the model, choices can be 'env',
% 'henv', 'ienv', 'penv', 'senv' and 'xenv'.
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
% *u*: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace for method 'env', 'henv', 'ienv', 'penv',
% ' senv', and 'xenv'.

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1 : 6);
%         modelType = 'env';
%         u = modelselectaic(X, Y, modelType)
% 
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X.X1 = fiberpaper(:, 7);
%         X.X2 = fiberpaper(:, 5 : 6);
%         modelType = 'penv';
%         u = modelselectaic(X, Y, modelType)

function u = modelselectaic(X, Y, modelType, Opts)

% Verify and initialize the parameters
%
if nargin < 3
    error('Inputs: X, Y and modelType should be specified!');
elseif nargin == 3
    Opts = [];
end

switch(modelType)
    case 'env'
        u = aic_env(X, Y, Opts);
    case 'henv'
        u = aic_henv(X, Y, Opts);
    case 'ienv'
        u = aic_ienv(X, Y, Opts);
    case 'penv'
        u = aic_penv(X, Y, Opts);
    case 'senv'
        u = aic_senv(X, Y, Opts);
    case 'xenv'
        u = aic_xenv(X, Y, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end