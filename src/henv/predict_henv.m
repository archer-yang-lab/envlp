%% predict_henv
% Perform estimation or prediction under the heteroscedastic envelope model.

%% Syntax
% PredictOutput=predict_henv(ModelOutput,Xnew,infType)
% PredictOutput=predict_henv(ModelOutput,Xnew,infType,Opts)
%
% Input
%
% ModelOutput: A list containing the maximum likelihood estimators and other
% statistics inherted from ienv.
% 
% Xnew: A group indicator.  It must be a column vector, whose transpose is
% the same as one of the group indictors from the original data. 
% 
% infType: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
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
% This function evaluates the inner envelope model at new value Xnew.  It can
% perform estimation: find the group mean for the group indicated by Xnew, or prediction:
% predict Y for the group indicated by Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
% load waterstrider.mat
% u=lrt_henv(X,Y,0.01);
% ModelOutput=henv(X,Y,u);
% ModelOutput.groupInd
% ModelOutput.mug
% Xnew=X(1,:)'
% PredictOutput=predict_henv(ModelOutput,Xnew,'estimation')
% PredictOutput.value  
% PredictOutput=predict_henv(ModelOutput,Xnew,'prediction')


function PredictOutput=predict_henv(ModelOutput,Xnew,infType,Opts)

if (nargin < 3)
    error('Inputs: ModelOutput,Xnew and infType should be specified!');
elseif (nargin==3)
    Opts=[];
end

if (~strcmp(infType,'estimation'))&&(~strcmp(infType,'prediction'))
    error('Inference type can only be estimation or prediction.');
end

[tmp iX iG]=intersect(Xnew',ModelOutput.groupInd,'rows');
if size(tmp,1)==0
    error('Xnew should be the same with one of the group indicators.')
end

[r u]=size(ModelOutput.Gamma);
ng=ModelOutput.ng;
n=sum(ng);

if u == 0
    
    if (strcmp(infType,'estimation'))
        
        PredictOutput.value=ModelOutput.mu;
        PredictOutput.covMatrix=ModelOutput.covMatrix(1:r)/n;
        PredictOutput.SE=sqrt(diag(PredictOutput.covMatrix));
        
    elseif (strcmp(infType,'prediction'))
        
        PredictOutput.value=ModelOutput.mug(:,iG);
        PredictOutput.covMatrix=ModelOutput.covMatrix(1:r)/n+ModelOutput.Sigma(:,:,iG);
        PredictOutput.SE=sqrt(diag(PredictOutput.covMatrix));
    
    end
       
else
    
    if (strcmp(infType,'estimation'))
        
        PredictOutput.value=ModelOutput.mug(:,iG);
        PredictOutput.covMatrix=(ModelOutput.covMatrix(1:r,1:r)+ModelOutput.covMatrix(iG*r+1:(iG+1)*r,iG*r+1:(iG+1)*r)+ModelOutput.covMatrix(1:r,iG*r+1:(iG+1)*r)+ModelOutput.covMatrix(iG*r+1:(iG+1)*r,1:r))/ng(iG);
        PredictOutput.SE=sqrt(diag(PredictOutput.covMatrix));
            
    elseif (strcmp(infType,'prediction'))
        
        PredictOutput.value=ModelOutput.mug(:,iG);
        PredictOutput.covMatrix=(ModelOutput.covMatrix(1:r,1:r)+ModelOutput.covMatrix(iG*r+1:(iG+1)*r,iG*r+1:(iG+1)*r)+ModelOutput.covMatrix(1:r,iG*r+1:(iG+1)*r)+ModelOutput.covMatrix(iG*r+1:(iG+1)*r,1:r))/ng(iG)+ModelOutput.Sigma(:,:,iG);
        PredictOutput.SE=sqrt(diag(PredictOutput.covMatrix));
    
        
    end
    
end