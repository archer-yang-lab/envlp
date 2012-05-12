%% predict_penv
% Perform estimation or prediction under the partial envelope model.

%% Syntax
% PredictOutput = predict_penv(ModelOutput, Xnew, infType)
%
%% Input
%
% *ModelOutput*: A list containing the maximum likelihood estimators and other
% statistics inherted from penv.
% 
% *Xnew*: A list containing the value of X1 and X2 with which to estimate or
% predict Y. 
% 
%  * Xnew.X1: A p1 by 1 vector containing the value of X1.
%  * Xnew.X2: A p2 by 1 vector containing the value of X2.
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
% This function evaluates the envelope model at new value Xnew.  It can
% perform estimation: find the fitted value when X = Xnew, or prediction:
% predict Y when X = Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         Xtemp = fiberpaper(:, 5 : 7);
%         X.X1 = Xtemp(:, 3);
%         X.X2 = Xtemp(:, 1 : 2);
%         alpha = 0.01;
%         u = lrt_penv(X, Y, alpha)
%         ModelOutput = penv(X, Y, u)
%         Xnew.X1 = X.X1(1, :)';
%         Xnew.X2 = X.X2(1, :)';
%         PredictOutput = predict_penv(ModelOutput, Xnew, 'estimation')
%         PredictOutput.value  % Compare the fitted value with the data
%         Y(1, :)'
%         PredictOutput = predict_penv(ModelOutput, Xnew, 'prediction')

function PredictOutput = predict_penv(ModelOutput, Xnew, infType)

if nargin < 3
    error('Inputs: ModelOutput, Xnew and infType should be specified!');
end

if ~strcmp(infType, 'estimation') && ~strcmp(infType, 'prediction')
    error('Inference type can only be estimation or prediction.');
end

[r p1] = size(ModelOutput.beta1);
p2 = size(ModelOutput.beta2, 2);

[s1 s2] = size(Xnew.X1);
if s1 ~= p1 || s2 ~= 1
    error('Xnew.X1 must be a p1 by 1 vector');
end

[s1 s2] = size(Xnew.X2);
if s1 ~= p2 || s2 ~= 1
    error('Xnew.X2 must be a p2 by 1 vector');
end

n = ModelOutput.n;
u = size(ModelOutput.Gamma, 2);

if u == 0
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.alpha;
        PredictOutput.covMatrix = ModelOutput.Sigma / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.alpha;
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.Sigma;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end
    
else
    
    if strcmp(infType, 'estimation')
        
        PredictOutput.value = ModelOutput.alpha + ModelOutput.beta1 * Xnew.X1 ...
            + ModelOutput.beta2 * Xnew.X2;
        X = [Xnew.X2' Xnew.X1']';
        PredictOutput.covMatrix = ModelOutput.Sigma / n ...
            + kron(X', eye(r)) * ModelOutput.covMatrix * kron(X, eye(r)) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    elseif strcmp(infType, 'prediction')
        
        PredictOutput.value = ModelOutput.alpha + ModelOutput.beta1 * Xnew.X1 ...
            + ModelOutput.beta2 * Xnew.X2;
        X = [Xnew.X2' Xnew.X1']';
        PredictOutput.covMatrix = (1 + 1 / n) * ModelOutput.Sigma ...
            + kron(X', eye(r)) * ModelOutput.covMatrix * kron(X, eye(r)) / n;
        PredictOutput.SE = sqrt(diag(PredictOutput.covMatrix));
        
    end

end