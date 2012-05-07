%% bstrp_senv
% Compute bootstrap standard error for the scaled envelope model. 

%% Syntax
% bootse=bstrp_senv(X,Y,u,B)
% bootse=bstrp_senv(X,Y,u,B,Opts)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% * u: Dimension of the envelope subspace.  A positive integer between 0 and
% r.
% * B: Number of boostrap samples.  A positive integer.
% * Opts: The optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used.
%
% Output
%
% * bootse: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the scaled envelope model by bootstrapping the residuals. 

%% Example
%
% load('T9-12.txt')
% Y=T9_12(:,4:7);
% X=T9_12(:,1:3);
% 
% u=bic_ienv(X,Y)
% B=20;
% bootse=bstrp_senv(X,Y,u,B)

function bootse=bstrp_senv(X,Y,u,B,Opts)

if (nargin < 4)
    error('Inputs: X, Y, u and B should be specified!');
elseif (nargin==4)
    Opts=[];
end

[n r]=size(Y);
p=size(X,2);

ModelOutput=senv(X,Y,u,Opts);

Yfit=ones(n,1)*ModelOutput.alpha'+X*ModelOutput.beta';
resi=Y-Yfit;

bootBeta=zeros(B,r*p);

for i=1:B
    
    bootresi=resi(randsample(1:n,n,true),:);
    Yboot=Yfit+bootresi;
    temp=senv(X,Yboot,u,Opts);
    bootBeta(i,:)=reshape(temp.beta,1,r*p);
    
end

bootse=reshape(sqrt(diag(cov(bootBeta,1))),r,p);