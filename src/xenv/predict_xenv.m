%% predict_xenv
% Perform estimation or prediction under the envelope model for the
% reduction on X.

%% Syntax
%         PredictOutput = predict_xenv(ModelOutput, Xnew, infType)
%
%% Input
%
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherted from xenv.
% 
% *Xnew*: The value of X with which to estimate or predict Y.  A p by 1
% vector.
% 
% *infType*: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
%
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
% This function evaluates the envelope model for the reduction on X at new value Xnew.  It can
% perform estimation: find the fitted value when X = Xnew, or prediction:
% predict Y when X = Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%         load wheatprotein.txt
%         X = wheatprotein(:, 1 : 6);
%         Y = wheatprotein(:, 7);
%         u = bic_xenv(X, Y)
%         ModelOutput = xenv(X, Y, u);
%         Xnew = X(1, :)';
%         PredictOutput = predict_xenv(ModelOutput, Xnew, 'estimation')
%         PredictOutput.value  % Compare the fitted value with the data
%         Y(1, :)'
%         PredictOutput = predict_xenv(ModelOutput, Xnew, 'prediction')

function PredictOutput = predict_xenv(ModelOutput, Xnew, infType)

if nargin < 3
    error('Inputs: ModelOutput, Xnew and infType should be specified!');
end

if ~strcmp(infType, 'estimation') && ~strcmp(infType, 'prediction')
    error('Inference type can only be estimation or prediction.');
end

[p r] = size(ModelOutput.beta);
[s1 s2] = size(Xnew);

if s1 ~= p || s2 ~= 1
    error('Xnew must be a p by 1 vector');
end

n = ModelOutput.n;
u = size(ModelOutput.Gamma, 2);

if u == 0
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = ModelOutput.sigYcX / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.sigYcX;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
else
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.mu + ModelOutput.beta' * Xnew;
        PredictOutput.covMatrix = ModelOutput.sigYcX / n ...
            + kron(eye(r), Xnew') * ModelOutput.covMatrix * kron(eye(r), Xnew) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.mu + ModelOutput.beta' * Xnew;
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.sigYcX ...
            + kron(eye(r), Xnew') * ModelOutput.covMatrix * kron(eye(r), Xnew) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
end