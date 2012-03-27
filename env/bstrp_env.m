%% bstrp_env
% Compute bootstrap standard error for the envelope model. 

%% Usage
% bootse=bstrp_env(X,Y,B,u)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% * B: Number of boostrap samples.  A positive integer.
% * u: Dimension of the envelope subspace.  A positive integer between 0 and
% r.
%
% Output
%
% * bootse: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the envelope model by bootstrapping the residuals. 

function bootse=bstrp_env(X,Y,B,u)


[n r]=size(Y);
p=size(X,2);

stat=env(X,Y,u);

Yfit=ones(n,1)*stat.alpha'+X*stat.beta';
resi=Y-Yfit;

bootBeta=zeros(B,r*p);

for i=1:B
    
    bootresi=resi(randsample(1:n,n,true),:);
    Yboot=Yfit+bootresi;
    temp=env(X,Yboot,u);
    bootBeta(i,:)=reshape(temp.beta,1,r*p);
    
end

bootse=reshape(sqrt(diag(cov(bootBeta,1))),r,p);