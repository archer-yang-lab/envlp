%% predict_ienv
% Perform estimation or prediction under the inner envelope model.

%% Syntax
% predOutput=predict_ienv(ModelOutput,Xnew,infType)
% predOutput=predict_ienv(ModelOutput,Xnew,infType,opts)
%
% Input
%
% ModelOutput: A list containing the maximum likelihood estimators and other
% statistics inherted from ienv.
% 
% Xnew: The value of X with which to estimate or predict Y.  A p by 1
% vector.
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
% perform estimation: find the fitted value when X=Xnew, or prediction:
% predict Y when X=Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
% load irisf.mat
% d=bic_ienv(X,Y);
% ModelOutput=ienv(X,Y,d);
% Xnew=X(1,:)';
% predOutput=predict_ienv(ModelOutput,Xnew,'estimation')
% predOutput.value  % Compare the fitted value with the data
% Y(1,:)'
% predOutput=predict_ienv(ModelOutput,Xnew,'prediction')


function predOutput=predict_ienv(ModelOutput,Xnew,infType,opts)

if (nargin < 3)
    error('Inputs: ModelOutput,Xnew and infType should be specified!');
elseif (nargin==3)
    opts=[];
end

if (~strcmp(infType,'estimation'))&&(~strcmp(infType,'prediction'))
    error('Inference type can only be estimation or prediction.');
end

[r,p]=size(ModelOutput.beta);
[s1 s2]=size(Xnew);

if s1~=p ||s2~=1
    error('Xnew must be a p by 1 vector');
end

n=ModelOutput.n;

if (strcmp(infType,'estimation'))
    
    predOutput.value=ModelOutput.alpha+ModelOutput.beta*Xnew;
    predOutput.covMatrix=ModelOutput.Sigma/n+kron(Xnew',eye(r))*ModelOutput.covMatrix*kron(Xnew,eye(r))/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
elseif (strcmp(infType,'prediction'))
    
    predOutput.value=ModelOutput.alpha+ModelOutput.beta*Xnew;
    predOutput.covMatrix=(1+1/n)*ModelOutput.Sigma+kron(Xnew',eye(r))*ModelOutput.covMatrix*kron(Xnew,eye(r))/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));    
    
end
