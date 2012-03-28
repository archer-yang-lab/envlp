%% lrt_penv
% Select the dimension of the partial envelope subspace using likelihood 
% ratio testing.

%% Usage
% u=lrt_penv(X1,X2,Y,alpha)
%
% Input
%
% * X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.
% * alpha: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
%
% Output
%
% * u: Dimension of the partial envelope subspace. An integer between 0 and r.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the partial envelope subspace, with prespecified 
% significance level $$\alpha$.  

%% Example
% load T7-7.dat
% Y=T7_7(:,1:4);
% X=T7_7(:,5:7);
% X1=X(:,3);
% X2=X(:,1:2);
% alpha=0.01;
% u=lrt_penv(X1,X2,Y,alpha)

function u=lrt_penv(X1,X2,Y,alpha)

[n r]=size(Y);

stat0=penv(X1,X2,Y,r);


for i=0:r-1
%     i
        stat=penv(X1,X2,Y,i);
        chisq = -2*(stat.l-stat0.l);
        df=stat0.np-stat.np;
        
        if (chi2cdf(chisq,df) < (1-alpha))
            u=i;
            break;
        end
end

if (i== r-1) && chi2cdf(chisq,df) > (1-alpha)
    u=r;
end