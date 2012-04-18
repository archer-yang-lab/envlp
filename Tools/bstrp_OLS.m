%% bstrp_OLS
% Compute bootstrap standard error for ordinary least squares. 

%% Usage
% bootse=bstrp_OLS(X,Y,B)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% * B: Number of boostrap samples.  A positive integer.
%
% Output
%
% * bootse: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in ordinary least squares by bootstrapping the residuals.  

function bootse=bstrp_OLS(X,Y,B)

dataParameter=make_parameter(X,Y);
n=dataParameter.n;
p=dataParameter.p;
r=dataParameter.r;
mY=dataParameter.mY;
XC=dataParameter.XC;
betaOLS=dataParameter.betaOLS;

Yfit=ones(n,1)*mY'+XC*betaOLS';
resi=Y-Yfit;

bootBeta=zeros(B,r*p);

for i=1:B
    
    bootresi=resi(randsample(1:n,n,true),:);
    Yboot=Yfit+bootresi;
    [betaBoot sigBoot]=fit_OLS(X,Yboot);
    bootBeta(i,:)=reshape(betaBoot,1,r*p);
    
end

bootse=reshape(sqrt(diag(cov(bootBeta,1))),r,p);