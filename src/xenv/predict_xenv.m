%% predict_xenv
% Perform estimation or prediction under the envelope model for the
% reduction on X.

%% Usage
% predOutput=predict_xenv(stat,Xnew,infType)
%
% Input
%
% * stat: A list containing the maximum likelihood estimators and other
% statistics inherted from xenv.
% * Xnew: The value of X with which to estimate or predict Y.  A p by 1
% vector.
% * infType: A string of characters indicting the inference type,
% the choices can be 'estimation' or 'prediction'.
%
% Output
%
% predOutput: A list containing the results of the inference.
%
% * predOutput.value: The fitted value or the prediction value evaluated at
% Xnew. An r by 1 vector.
% * predOutput.covMatrix: The covariance matrix of predOutput.value. An r by r
% matrix.
% * predOutput.SE: The standard error of elements in predOutput.value.

%% Description
% This function evaluates the envelope model for the reduction on X at new value Xnew.  It can
% perform estimation: find the fitted value when X=Xnew, or prediction:
% predict Y when X=Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
% load wheatprotein.txt
% X=wheatprotein(:,1:6);
% Y=wheatprotein(:,7);
% u=bic_xenv(X,Y)
% stat=xenv(X,Y,u);
% Xnew=X(1,:)';
% predOutput=predict_xenv(stat,Xnew,'estimation')
% predOutput.value  % Compare the fitted value with the data
% Y(1,:)'
% predOutput=predict_senv(stat,Xnew,'prediction')

function predOutput=predict_xenv(stat,Xnew,infType)

if (nargin < 3)
    error('Inputs: stat,Xnew and infType should be specified!');
end

if (~strcmp(infType,'estimation'))&&(~strcmp(infType,'prediction'))
    error('Inference type can only be estimation or prediction.');
end

[p,r]=size(stat.beta);
[s1 s2]=size(Xnew);

if s1~=p ||s2~=1
    error('Xnew must be a p by 1 vector');
end

n=stat.n;

if (strcmp(infType,'estimation'))
    
    predOutput.value=stat.mu+stat.beta'*Xnew;
    predOutput.covMatrix=stat.sigYcX/n+kron(eye(r),Xnew')*stat.covMatrix*kron(eye(r),Xnew)/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
elseif (strcmp(infType,'prediction'))
    
    predOutput.value=stat.mu+stat.beta'*Xnew;
    predOutput.covMatrix=(1+1/n)*stat.sigYcX+kron(eye(r),Xnew')*stat.covMatrix*kron(eye(r),Xnew)/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
end
