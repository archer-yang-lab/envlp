%% predict_penv
% Perform estimation or prediction under the partial envelope model.

%% Syntax
% predOutput=predict_penv(ModelOutput,Xnew,infType)
% predOutput=predict_penv(ModelOutput,Xnew,infType,opts)
%
% Input
%
% ModelOutput: A list containing the maximum likelihood estimators and other
% statistics inherted from penv.
% 
% Xnew: A list containing the value of X1 and X2 with which to estimate or
% predict Y. 
% 
%  * Xnew.X1: A p1 by 1 vector containing the value of X1.
%  * Xnew.X2: A p2 by 1 vector containing the value of X2.
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
% This function evaluates the envelope model at new value Xnew.  It can
% perform estimation: find the fitted value when X=Xnew, or prediction:
% predict Y when X=Xnew.  The covariance matrix and the standard errors are
% also provided.

%% Example
%
% load T7-7.dat
% Y = T7_7(:,1:4);
% X = T7_7(:,5:7);
% X1 = X(:,3);
% X2 = X(:,1:2);
% alpha = 0.01;
% u = lrt_penv(X1,X2,Y,alpha)
% ModelOutput = penv(X1,X2,Y,u)
% Xnew.X1=X1(1,:)';
% Xnew.X2=X2(1,:)';
% predOutput=predict_penv(ModelOutput,Xnew,'estimation')
% predOutput.value  % Compare the fitted value with the data
% Y(1,:)'
% predOutput=predict_penv(ModelOutput,Xnew,'prediction')

function predOutput=predict_penv(ModelOutput,Xnew,infType,opts)

if (nargin < 3)
    error('Inputs: ModelOutput,Xnew and infType should be specified!');
elseif (nargin==3)
    opts=[];
end

if (~strcmp(infType,'estimation'))&&(~strcmp(infType,'prediction'))
    error('Inference type can only be estimation or prediction.');
end

[r,p1]=size(ModelOutput.beta1);
p2=size(ModelOutput.beta2,2);

[s1 s2]=size(Xnew.X1);
if s1~=p1 ||s2~=1
    error('Xnew.X1 must be a p1 by 1 vector');
end

[s1 s2]=size(Xnew.X2);
if s1~=p2 ||s2~=1
    error('Xnew.X2 must be a p2 by 1 vector');
end

n=ModelOutput.n;

if (strcmp(infType,'estimation'))
    
    predOutput.value=ModelOutput.alpha+ModelOutput.beta1*Xnew.X1+ModelOutput.beta2*Xnew.X2;
    X=[Xnew.X2' Xnew.X1']';
    predOutput.covMatrix=ModelOutput.Sigma/n+kron(X',eye(r))*ModelOutput.covMatrix*kron(X,eye(r))/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
elseif (strcmp(infType,'prediction'))
    
    predOutput.value=ModelOutput.alpha+ModelOutput.beta1*Xnew.X1+ModelOutput.beta2*Xnew.X2;
    X=[Xnew.X2' Xnew.X1']';
    predOutput.covMatrix=(1+1/n)*ModelOutput.Sigma+kron(X',eye(r))*ModelOutput.covMatrix*kron(X,eye(r))/n;
    predOutput.SE=sqrt(diag(predOutput.covMatrix));
    
end
