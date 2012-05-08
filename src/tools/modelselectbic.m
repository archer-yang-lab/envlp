%% modelselectbic
%
%% Syntax
% u = modelselectbic(X, Y, modelType)
% u = modelselectbic(X, Y, modelType, Opts)
% 
% Input
%
% X: Predictors.   The predictors can be univariate or multivariate, 
% discrete or continuous.  For model type for method 'env', 'henv', 'ienv',
% ' senv', and 'xenv'. X is an n by p matrix, p is the number of
% predictors. For model type 'penv', X is  A list containing the value of X1 and X2.
% 
% * X.X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
%
% Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% Opts: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out output, logical 0 or 1. Default value:
% 0.
%
% Output
%
% u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the envelope subspace for method 'env', 'henv', 'ienv', 'penv',
% ' senv', and 'xenv'.

%% Example
% load wheatprotein.txt
% X = wheatprotein(:, 8);
% Y = wheatprotein(:, 1 : 6);
% modelType = 'env';
% u = modelselectbic(X, Y, modelType)
% load T7-7.dat
% Y = T7_7(:, 1 : 4);
% Xtemp = T7_7(:, 5 : 7);
% X.X1 = Xtemp(:, 3);
% X.X2 = Xtemp(:, 1 : 2);
% modelType = 'penv';
% u = modelselectbic(X, Y, modelType)

function u = modelselectbic(X, Y, modelType, Opts)

% Verify and initialize the parameters
%
if nargin < 3
    error('Inputs: X, Y and modelType should be specified!');
elseif nargin == 3
    Opts = [];
end

switch(modelType)
    case 'env'
        u = bic_env(X, Y, Opts);
    case 'henv'
        u = bic_henv(X, Y, Opts);
    case 'ienv'
        u = bic_ienv(X, Y, Opts);
    case 'penv'
        u = bic_penv(X, Y, Opts);
    case 'senv'
        u = bic_senv(X, Y, Opts);
    case 'xenv'
        u = bic_xenv(X, Y, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end