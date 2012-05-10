%% bootstrapse
% Perform boostrap to estimate actual standard erros for models in the envelope family.
%
%% Syntax
%         bootse = bootstrapse(X, Y, u, B, modelType)
%         bootse = bootstrapse(X, Y, u, B, modelType, Opts)
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
% * X.X1 (only for 'penv'): Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2 (only for 'penv'): Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% *u*: Dimension of the envelope subspace. The legitimate range of u depends
% on the model specified. 
% 
% *B*: Number of boostrap samples.  A positive integer.
% 
% *modelType*: A string characters indicting the model, choices can be 'env',
% 'henv', 'ienv', 'penv', 'senv' and 'xenv'.
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
% bootse: For 'env', 'henv', 'ienv', 'senv' and 'xenv', an r by p matrix 
% containing the standard errors for elements in $$\beta$ computed by
% bootstrap.  For 'penv', an r by p1 matrix containing the standard errors 
% for $$\beta_1$ computed by bootstrap.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients or for partial envelope model, the main regression 
% coefficients in the specified model by bootstrapping the residuals.

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1:6);
%         alpha = 0.01;
%         u = lrt_env(X, Y, alpha)
%         B = 100;
%         modelType = 'env';
%         bootse = bootstrapse(X, Y, u, B, modelType)
% 
%         load T7-7.dat
%         Y = T7_7(:, 1 : 4);
%         Xtemp = T7_7(:, 5 : 7);
%         X.X1 = Xtemp(:, 3);
%         X.X2 = Xtemp(:, 1 : 2);
%         alpha = 0.01;
%         u = lrt_penv(X, Y, alpha);
%         B = 100;
%         modelType = 'penv';
%         bootse = bootstrapse(X, Y, u, B, modelType)
% 

function bootse = bootstrapse(X, Y, u, B, modelType, Opts)

% Verify and initialize the parameters
%
if nargin < 5
    error('Inputs: X, Y, u, B and modelType should be specified!');
elseif nargin == 5
    Opts = [];
end

switch(modelType)
    case 'env'
        bootse = bstrp_env(X, Y, u, B, Opts);
    case 'henv'
        bootse = bstrp_henv(X, Y, u, B, Opts);
    case 'ienv'
        bootse = bstrp_ienv(X, Y, u, B, Opts);
    case 'penv'
        bootse = bstrp_penv(X, Y, u, B, Opts);
    case 'senv'
        bootse = bstrp_senv(X, Y, u, B, Opts);
    case 'xenv'
        bootse = bstrp_xenv(X, Y, u, B, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end