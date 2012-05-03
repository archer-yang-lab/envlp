%% bstrp_ienv
% Compute bootstrap standard error for the inner envelope model. 

%% Usage
% bootse=bstrp_ienv(X,Y,u,B,opts)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% * u: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.
% * B: Number of boostrap samples.  A positive integer.
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
% bootse=bstrp_ienv(X,Y,u,B)

function bootse=bstrp_ienv(X,Y,u,B,opts)

if (nargin < 4)
    error('Inputs: X, Y, B and u should be specified!');
elseif (nargin==4)
    opts=[];
end

[n r]=size(Y);
p=size(X,2);

stat=ienv(X,Y,u,opts);

Yfit=ones(n,1)*stat.alpha'+X*stat.beta';
resi=Y-Yfit;

bootBeta=zeros(B,r*p);

for i=1:B
    
    bootresi=resi(randsample(1:n,n,true),:);
    Yboot=Yfit+bootresi;
    temp=ienv(X,Yboot,u,opts);
    bootBeta(i,:)=reshape(temp.beta,1,r*p);
    
end

bootse=reshape(sqrt(diag(cov(bootBeta,1))),r,p);