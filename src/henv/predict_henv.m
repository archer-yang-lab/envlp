%% predict_henv
% Perform estimation or prediction under the heteroscedastic envelope model.

%% Syntax
%         PredictOutput = predict_henv(ModelOutput, Xnew, infType)
%
%% Input
%
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherited from henv.
% 
% *Xnew*: A group indicator.  It must be a column vector, whose transpose is
% the same as one of the group indictors from the original data. 
% 
% *infType*: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
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
% This function evaluates the inner envelope model at new value Xnew.  It can
% perform estimation: find the group mean for the group indicated by Xnew, or prediction:
% predict Y for the group indicated by Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
%         load waterstrider.mat
%         u = lrt_henv(X, Y, 0.01);
%         ModelOutput = henv(X, Y, u);
%         ModelOutput.groupInd
%         ModelOutput.mug
%         Xnew = X(1, :)'
%         PredictOutput = predict_henv(ModelOutput, Xnew, 'estimation')
%         PredictOutput.value %This is the 3rd group mean
%         PredictOutput.SE
%         PredictOutput = predict_henv(ModelOutput, Xnew, 'prediction')
%         PredictOutput.SE


function PredictOutput = predict_henv(ModelOutput, Xnew, infType)

if nargin < 3
    error('Inputs: ModelOutput, Xnew and infType should be specified!');
end

if ~strcmp(infType, 'estimation') && ~strcmp(infType, 'prediction')
    error('Inference type can only be estimation or prediction.');
end

[tmp, ~, iG] = intersect(Xnew', ModelOutput.groupInd, 'rows');
if size(tmp, 1)==0
    error('Xnew should be the same with one of the group indicators.')
end

[r, u] = size(ModelOutput.Gamma);
ng = ModelOutput.ng;
n = sum(ng);

if u == 0
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.mu;
        PredictOutput.covMatrix = ModelOutput.covMatrix(1:r) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.mug(:, iG);
        PredictOutput.covMatrix = ModelOutput.covMatrix(1:r) / n + ModelOutput.Sigma(:, :, iG);
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
    
    end
       
else
    
    temp1 = iG * r + 1;
    temp2 = (iG + 1) * r;
    
    if strcmp(infType, 'estimation')
                
        PredictOutput.value = ModelOutput.mug(:, iG);
        PredictOutput.covMatrix = (ModelOutput.covMatrix(1 : r, 1 : r) ...
            + ModelOutput.covMatrix(temp1 : temp2, temp1 : temp2) ...
            + ModelOutput.covMatrix(1 : r, temp1 : temp2) ...
            + ModelOutput.covMatrix(temp1 : temp2, 1 : r)) / ng(iG);
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
            
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.mug(:, iG);
        PredictOutput.covMatrix = (ModelOutput.covMatrix(1 : r, 1 : r) ...
            + ModelOutput.covMatrix(temp1 : temp2, temp1 : temp2) ...
            + ModelOutput.covMatrix(1 : r, temp1 : temp2) ...
            + ModelOutput.covMatrix(temp1 : temp2, 1 : r)) / ng(iG) ...
            + ModelOutput.Sigma(:, :, iG);
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
           
    end
    
end