%% mfoldcv
% Select the dimension for the envelope family using m-fold cross
% validation.

%% Syntax
%         u = mfoldcv(X, Y, m, modelType)
%         u = mfoldcv(X, Y, m, modelType, Opts)
% 
%% Input
%
% *X*: Predictors.  The predictors can be univariate or multivariate, 
% discrete or continuous.  
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% 
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *modelType*: A string characters indicting the model, choices can be
% 'envseq' or 'xenvpls'.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out the number of bootstrap samples, 
% logical 0 or 1. Default value: 0.

%% Output
%
%  *u*: The dimension of the envelope subspace selected by m-fold cross
%  validation.

%% Description
% This function implements m-fold cross validation to select the dimension
% of the envelope space, based on prediction performance.  For each u, the
% data is partitioned into m parts, each part is in turn used for testing 
% for the prediction performance while the rest m-1 parts are used for 
% training.  The dimension is select as the one that minimizes the average 
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
    case 'envseq'
        u = mfoldcv_envseq(X, Y, m, Opts);
    case 'xenvpls'
        u = mfoldcv_xenvpls(X, Y, m, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end