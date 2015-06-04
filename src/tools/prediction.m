%% prediction
% Perform estimation or prediction for models in the envelope family.
%
%% Syntax
%         PredictOutput = prediction(ModelOutput, Xnew, infType, modelType)
% 
%% Input
%
% *ModelOutput*: A list containing the model outputs from fitting the models.
%
% *Xnew*: The value of X with which to estimate or predict Y.  
% 
% For 'env',
% 'henv', 'ienv', 'senv', 'sxenv' and 'xenv', it is a p by 1 vector.  
% 
% For 'penv', it is a list containing the value of X1 and X2.
% 
%  * Xnew.X1 (only for 'penv'): A p1 by 1 vector containing the value of X1.
%  * Xnew.X2 (only for 'penv'): A p2 by 1 vector containing the value of X2.
% 
% *infType*: A string of characters indicating the inference type,
% the choices can be 'estimation' or 'prediction'.
% 
% *modelType*: A string of characters indicating the model, choices can be 'env',
% 'henv', 'ienv', 'penv', 'senv', 'sxenv' and 'xenv'.
%
%% Output
%
% *PredictOutput*: A list containing the results of the inference.
%
% * PredictOutput.value: The fitted value or the prediction value evaluated at
% Xnew. An r by 1 vector.
% * PredictOutput.covMatrix: The covariance matrix of PredictOutput.value. An r by r
% matrix.
% * PredictOutput.SE: The standard error of elements in PredictOutput.value. An r
% by 1 vector. 

%% Description
% This function evaluates the user-specified model, could be 'env', 'henv', 'ienv', 'penv', 
% 'senv', 'sxenv' or 'xenv', at new value Xnew.  It can
% perform estimation: find the fitted value when X = Xnew, or prediction:
% predict Y when X = Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 8);
%         Y = wheatprotein(:, 1:6);
%         modelType = 'env';
%         u =  modelselectbic(X, Y, modelType);
%         ModelOutput = env(X, Y, u);
%         Xnew = X(2, :)';
%         PredictOutput = prediction(ModelOutput, Xnew, 'estimation', modelType)
%         [PredictOutput.value, Y(2, :)'] % Compare the fitted value with the observed value
% 
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X.X1 = fiberpaper(:, 7);
%         X.X2 = fiberpaper(:, 5 : 6);
%         modelType = 'penv';
%         u =  modelselectbic(X, Y, modelType);
%         ModelOutput = penv(X, Y, u);
%         Xnew.X1 = X.X1(1, :)';
%         Xnew.X2 = X.X2(1, :)';
%         PredictOutput = prediction(ModelOutput, Xnew, 'estimation', modelType)
%         PredictOutput.SE
% 

function PredictOutput = prediction(ModelOutput, Xnew, infType, modelType)

% Verify and initialize the parameters
%
if nargin < 4
    error('Inputs: ModelOutput, Xnew, infType and modelType should be specified!');
end

switch(modelType)
    case 'env'
        PredictOutput = predict_env(ModelOutput, Xnew, infType);
    case 'henv'
        PredictOutput = predict_henv(ModelOutput, Xnew, infType);
    case 'ienv'
        PredictOutput = predict_ienv(ModelOutput, Xnew, infType);
    case 'penv'
        PredictOutput = predict_penv(ModelOutput, Xnew, infType);
    case 'senv'
        PredictOutput = predict_senv(ModelOutput, Xnew, infType);
    case 'sxenv'
        PredictOutput = predict_sxenv(ModelOutput, Xnew, infType);
    case 'xenv'
        PredictOutput = predict_xenv(ModelOutput, Xnew, infType);
    otherwise
        fprintf('The value specified in modelType is not supported!');
end