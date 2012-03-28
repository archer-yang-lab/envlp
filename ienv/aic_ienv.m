%% aic_ienv
% Select the dimension of the inner envelope subspace using Akaike
% information criterion.

%% Usage
% u=aic_ienv(X,Y)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors and n 
% is the number of observations. The predictors can be univariate or 
% multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.
%
% Output
%
% * u: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the inner envelope subspace. 

%% Example
%
% load irisf.mat
% u=aic_ienv(X,Y)

function u=aic_ienv(X,Y)

[n r]=size(Y);
p=size(X,2);
    
stat=ienv(X,Y,r);
ic=-2*stat.l+2*stat.np;
u=r;


for i=0:p
%     i
        stat=ienv(X,Y,i);
        temp=-2*stat.l+2*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end