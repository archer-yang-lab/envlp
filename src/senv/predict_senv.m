%% predict_senv
% Perform estimation or prediction under the scaled envelope model.

%% Usage
% predOutput=predict_senv(stat,Xnew,infType)
%
% Input
%
% * stat: A list containing the maximum likelihood estimators and other
% statistics inherted from senv.
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
% This function evaluates the scaled envelope model at new value Xnew.  It can
% perform estimation: find the fitted value when X=Xnew, or prediction:
% predict Y when X=Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
% load('T9-12.txt')
% Y=T9_12(:,4:7);
% X=T9_12(:,1:3);
% u=bic_env(X,Y)
% stat=env(X,Y,u);
% Xnew=X(1,:)';
% predOutput=predict_senv(stat,Xnew,'estimation')
% predOutput.value  % Compare the fitted value with the data
% Y(1,:)'
% predOutput=predict_senv(stat,Xnew,'prediction')

function predOutput=predict_senv(stat,Xnew,infType)

if (nargin < 3)
    error('Inputs: stat,Xnew and infType should be specified!');
end

if (~strcmp(infType,'estimation'))&&(~strcmp(infType,'prediction'))
    error('Inference type can only be estimation or prediction.');
end

[r,p]=size(stat.beta);
[s1 s2]=size(Xnew);

if s1~=p ||s2~=1
    error('Xnew must be a p by 1 vector');
end

n=stat.n;

if (strcmp(infType,'estimation'))
    
    predOutput.value=stat.alpha+stat.beta*Xnew;
    predOutput.covMatrix=stat.Sigma/n+kron(Xnew',eye(r))*stat.covMatrix*kron(Xnew,eye(r))/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
elseif (strcmp(infType,'prediction'))
    
    predOutput.value=stat.alpha+stat.beta*Xnew;
    predOutput.covMatrix=(1+1/n)*stat.Sigma+kron(Xnew',eye(r))*stat.covMatrix*kron(Xnew,eye(r))/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
end
