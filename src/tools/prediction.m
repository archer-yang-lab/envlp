%% prediction
% Perform estimation or prediction for models in the envelope family.
%
%% Syntax
% PredictOutput = prediction(ModelOutput, Xnew, infType, modelType)
% PredictOutput = prediction(ModelOutput, Xnew, infType, modelType, Opts)
% 
% Input
%
% ModelOutput: A list containing the model outputs from fitting the models.
%
% Xnew: The value of X with which to estimate or predict Y.  
% 
% For 'env',
% 'henv', 'ienv', 'senv' and 'xenv', it is a p by 1 vector.  
% 
% For 'penv', it is a list containing the value of X1 and X2.
% 
%  * Xnew.X1 (only for 'penv'): A p1 by 1 vector containing the value of X1.
%  * Xnew.X2 (only for 'penv'): A p2 by 1 vector containing the value of X2.
% 
% infType: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
% 
% modelType: A string characters indicting the model, choices can be 'env',
% 'henv', 'ienv', 'penv', 'senv' and 'xenv'.
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
% PredictOutput: A list containing the results of the inference.
%
% * PredictOutput.value: The fitted value or the prediction value evaluated at
% Xnew. An r by 1 vector.
% * PredictOutput.covMatrix: The covariance matrix of PredictOutput.value. An r by r
% matrix.
% * PredictOutput.SE: The standard error of elements in PredictOutput.value. An r
% by 1 vector. 

%% Description
% This function evaluates the model, could be 'env', 'henv', 'ienv', 'penv', 
% 'senv' or 'xenv', at new value Xnew.  It can
% perform estimation: find the fitted value when X = Xnew, or prediction:
% predict Y when X = Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
% load wheatprotein.txt
% X = wheatprotein(:, 8);
% Y = wheatprotein(:, 1:6);
% modelType = 'env';
% u =  modelselectbic(X, Y, modelType)
% ModelOutput = env(X, Y, u);
% Xnew = X(2, :)';
% PredictOutput = predict_env(ModelOutput, Xnew, 'estimation')
% 
% load T7-7.dat
% Y = T7_7(:, 1 : 4);
% Xtemp = T7_7(:, 5 : 7);
% X.X1 = Xtemp(:, 3);
% X.X2 = Xtemp(:, 1 : 2);
% modelType = 'penv';
% u =  modelselectbic(X, Y, modelType)
% ModelOutput = penv(X, Y, u)
% Xnew.X1 = X.X1(1, :)';
% Xnew.X2 = X.X2(1, :)';
% PredictOutput = predict_penv(ModelOutput, Xnew, 'estimation')


function PredictOutput = prediction(ModelOutput, Xnew, infType, modelType, Opts)

% Verify and initialize the parameters
%
if nargin < 4
    error('Inputs: ModelOutput, Xnew, infType and modelType should be specified!');
elseif nargin == 4
    Opts = [];
end

switch(modelType)
    case 'env'
        u = predict_env(ModelOutput, Xnew, infType, Opts);
    case 'henv'
        u = predict_henv(ModelOutput, Xnew, infType, Opts);
    case 'ienv'
        u = predict_ienv(ModelOutput, Xnew, infType, Opts);
    case 'penv'
        u = predict_penv(ModelOutput, Xnew, infType, Opts);
    case 'senv'
        u = predict_senv(ModelOutput, Xnew, infType, Opts);
    case 'xenv'
        u = predict_xenv(ModelOutput, Xnew, infType, Opts);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end