%% predict_envmean
% Perform estimation of the multivariate mean or prediction for a new observation.

%% Syntax
%         PredictOutput = predict_envmean(ModelOutput, infType)
%
%% Input
%
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherited from envmean.
% 
% *infType*: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
% 
%% Output
%
% *PredictOutput*: A list containing the results of the inference.
%
% * PredictOutput.value: The estimated multivariate mean or the prediction
% value. A p dimensional column vector. 
% * PredictOutput.covMatrix: The covariance matrix of PredictOutput.value.
% A p by p matrix.
% * PredictOutput.SE: The standard error of elements in
% PredictOutput.value.  A p dimensional column vector. 

%% Description
% If the inference type is prediction, this function predicts a new
% observation and gives its covariance matrix and standard errors of its elements.  If the
% inference type is estimation, this function gives the estimation of the
% multivariate mean, its covariance matrix and standard errors of its elements.

%% Example
%
%         load Adopted
%         Y = Adopted(:, 1 : 6);
%         u = bic_envmean(Y);
%         ModelOutput = envmean(Y, u);
%         PredictOutput = predict_envmean(ModelOutput, 'prediction')
%         PredictOutput.value
%         PredictOutput.SE

function PredictOutput = predict_envmean(ModelOutput, infType)

if nargin < 2
    error('Inputs: ModelOutput and infType should be specified!');
end

if ~strcmp(infType, 'estimation') && ~strcmp(infType, 'prediction')
    error('Inference type can only be estimation or prediction.');
end


n = ModelOutput.n;
u = size(ModelOutput.Gamma, 2);
p = size(ModelOutput.Sigma, 1);

if u == 0
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = zeros(p, p);
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = ModelOutput.Sigma;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
else
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = ModelOutput.Sigma / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.Sigma;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
end