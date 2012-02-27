%% aic_env
% Select the dimension of the envelope subspace using Akaike information
% criterion.

%% Usage
% u=aic_env(Y,X)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% Output
%
% * u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace.  

function u=aic_env(Y,X)

[n r]=size(Y);
p=size(X,2);
    
stat=env(X,Y,r);
ic=-2*stat.l+2*stat.np;
u=r;


for i=0:r-1
%     i
        stat=env(X,Y,i);
        temp=-2*stat.l+2*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end
