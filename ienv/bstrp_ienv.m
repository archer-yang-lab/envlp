%% bstrp_ienv
% Compute bootstrap standard error for the inner envelope model. 

%% Usage
% bootse=bstrp_ienv(X,Y,B,u)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% * B: Number of boostrap samples.  A positive integer.
% * u: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.
%
% Output
%
% * bootse: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the inner envelope model by bootstrapping the residuals. 

%% Example
%
% load irisf.mat
% 
% u=bic_ienv(X,Y)
% B=100;
% bootse=bstrp_ienv(X,Y,B,u)

function bootse=bstrp_ienv(X,Y,B,u)


[n r]=size(Y);
p=size(X,2);

stat=ienv(X,Y,u);

Yfit=ones(n,1)*stat.alpha'+X*stat.beta';
resi=Y-Yfit;

bootBeta=zeros(B,r*p);

for i=1:B
    
    bootresi=resi(randsample(1:n,n,true),:);
    Yboot=Yfit+bootresi;
    temp=ienv(X,Yboot,u);
    bootBeta(i,:)=reshape(temp.beta,1,r*p);
    
end

bootse=reshape(sqrt(diag(cov(bootBeta,1))),r,p);