%% predict_henv
% Perform estimation or prediction under the heteroscedastic envelope model.

%% Syntax
% predOutput=predict_henv(stat,Xnew,infType)
% predOutput=predict_henv(stat,Xnew,infType,opts)
%
% Input
%
% stat: A list containing the maximum likelihood estimators and other
% statistics inherted from ienv.
% 
% Xnew: A group indicator.  It must be a column vector, whose transpose is
% the same as one of the group indictors from the original data. 
% 
% infType: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
%
% opts: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * opts.maxIter: Maximum number of iterations.  Default value: 300.
% * opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * opts.verbose: Flag for print out output, logical 0 or 1. Default value:
% 0.
% 
% Output
%
% predOutput: A list containing the results of the inference.
%
% * predOutput.value: The fitted value or the prediction value evaluated at
% Xnew. An r by 1 vector.
% * predOutput.covMatrix: The covariance matrix of predOutput.value. An r by r
% matrix.
% * predOutput.SE: The standard error of elements in predOutput.value. An r
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
% stat=henv(X,Y,u);
% stat.groupInd
% stat.mug
% Xnew=X(1,:)'
% predOutput=predict_henv(stat,Xnew,'estimation')
% predOutput.value  
% predOutput=predict_henv(stat,Xnew,'prediction')


function predOutput=predict_henv(stat,Xnew,infType,opts)

if (nargin < 3)
    error('Inputs: stat,Xnew and infType should be specified!');
elseif (nargin==3)
    opts=[];
end

if (~strcmp(infType,'estimation'))&&(~strcmp(infType,'prediction'))
    error('Inference type can only be estimation or prediction.');
end

[tmp iX iG]=intersect(Xnew',stat.groupInd,'rows');
if size(tmp,1)==0
    error('Xnew should be the same with one of the group indicators.')
end

[r u]=size(stat.Gamma);
ng=stat.ng;
n=sum(ng);

if u == 0
    
    if (strcmp(infType,'estimation'))
        
        predOutput.value=stat.mu;
        predOutput.covMatrix=stat.covMatrix(1:r)/n;
        predOutput.SE=sqrt(diag(predOutput.covMatrix));
        
    elseif (strcmp(infType,'prediction'))
        
        predOutput.value=stat.mug(:,iG);
        predOutput.covMatrix=stat.covMatrix(1:r)/n+stat.Sigma(:,:,iG);
        predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
    end
       
else
    
    if (strcmp(infType,'estimation'))
        
        predOutput.value=stat.mug(:,iG);
        predOutput.covMatrix=(stat.covMatrix(1:r,1:r)+stat.covMatrix(iG*r+1:(iG+1)*r,iG*r+1:(iG+1)*r)+stat.covMatrix(1:r,iG*r+1:(iG+1)*r)+stat.covMatrix(iG*r+1:(iG+1)*r,1:r))/ng(iG);
        predOutput.SE=sqrt(diag(predOutput.covMatrix));
            
    elseif (strcmp(infType,'prediction'))
        
        predOutput.value=stat.mug(:,iG);
        predOutput.covMatrix=(stat.covMatrix(1:r,1:r)+stat.covMatrix(iG*r+1:(iG+1)*r,iG*r+1:(iG+1)*r)+stat.covMatrix(1:r,iG*r+1:(iG+1)*r)+stat.covMatrix(iG*r+1:(iG+1)*r,1:r))/ng(iG)+stat.Sigma(:,:,iG);
        predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
        
    end
    
end