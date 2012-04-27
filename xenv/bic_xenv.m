%% bic_xenv
% Use Bayesian information criterion to select the dimension of the envelope
% subspace for the reduction on X.

%% Usage
% u=bic_xenv(X,Y)
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
% * u: Dimension of the envelope. An integer between 0 and p.

%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the envelope subspace for the reduction on X.

%% Example
%
% load wheatprotein.txt
% X=wheatprotein(:,1:6);
% Y=wheatprotein(:,7);
% u=bic_xenv(X,Y)

function u=bic_xenv(X,Y)

[n p]=size(X);
    
stat=xenv(X,Y,p);
ic=-2*stat.l+log(n)*stat.np;
u=p;


for i=0:p-1
%     i
        stat=xenv(X,Y,i);
        temp=-2*stat.l+2*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end
