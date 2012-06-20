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