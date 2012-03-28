%% bic_penv
% Select the dimension of the partial envelope subspace using Bayesian information
% criterion.

%% Usage
% u=bic_penv(X1,X2,Y)
%
% Input
%
% * X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% Output
%
% * u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the partial envelope subspace. 

%% Example
% load T7-7.dat
% Y=T7_7(:,1:4);
% X=T7_7(:,5:7);
% X1=X(:,3);
% X2=X(:,1:2);
% u=bic_penv(X1,X2,Y)

function u=bic_penv(X1,X2,Y)

[n r]=size(Y);
    
stat=penv(X1,X2,Y,r);
ic=-2*stat.l+log(n)*stat.np;
u=r;


for i=0:r-1
%     i
        stat=penv(X1,X2,Y,i);
        temp=-2*stat.l+log(n)*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end
