%% mfoldcv
% Select the dimension for the envelope family using m-fold cross
% validation.

%% Syntax
%         u = mfoldcv(X, Y, m, modelType)
%         u = mfoldcv(X, Y, m, modelType, Opts)
% 
%% Input
%
% *X*: Predictors.   The predictors can be univariate or multivariate, 
% discrete or continuous.  
% 
% For model type 'env', 'envseq', 'henv', 'ienv',
% ' senv', 'xenv' and 'xenvpls', X is an n by p matrix, p is the number of
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
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *modelType*: A string of characters indicting the model, choices can be
% 'env', 'envseq', 'henv', 'ienv', 'penv', 'senv', 'xenv' or 'xenvpls'.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag to print out dimension selection process, 
% logical 0 or 1. Default value: 0.
% * Opts.table: Flag to tabulate the results, which contains cross 
% validation error for each u.  Logical 0 or 1. Default value: 0.

%% Output
%
%  *u*: The dimension of the envelope subspace selected by m-fold cross
%  validation.

%% Description
% This function implements m-fold cross validation to select the dimension
% of the envelope space, based on prediction performance.  For each u, the
% data is partitioned into m parts, each part is in turn used for testing 
% for the prediction performance while the rest m-1 parts are used for 
% training.  The dimension is selected as the one that minimizes the average 
% prediction errors. If Y is multivariate, the identity inner product is 
% used for computing the prediction errors.

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1 : 6);
%         m = 5;
%         modelType = 'envseq';
%         u = mfoldcv(X, Y, m, modelType)


function u = mfoldcv(X, Y, m, modelType, Opts)

% Verify and initialize the parameters
%
if nargin < 4
    error('Inputs: X, Y, m and modelType should be specified!');
elseif nargin == 4
    Opts = [];
end

switch(modelType)
    case 'env'
        u = mfoldcv_env(X, Y, m, Opts);
    case 'envseq'
        u = mfoldcv_envseq(X, Y, m, Opts);
    case 'henv'
        u = mfoldcv_henv(X, Y, m, Opts);
    case 'ienv'
        u = mfoldcv_ienv(X, Y, m, Opts);
    case 'penv'
        u = mfoldcv_penv(X, Y, m, Opts);
    case 'senv'
        u = mfoldcv_senv(X, Y, m, Opts);
    case 'xenv'
        u = mfoldcv_xenv(X, Y, m, Opts);
    case 'xenvpls'
        u = mfoldcv_xenvpls(X, Y, m, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end