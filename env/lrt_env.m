%% lrt_env
% Select the dimension of the envelope subspace using likelihood ratio
% testing.

%% Usage
% u=lrt_env(X,Y,alpha)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% * alpha: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
%
% Output
%
% * u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the envelope subspace, with prespecified significance 
% level $$\alpha$.  

%% Example
% load wheatprotein.txt
% X=wheatprotein(:,8);
% Y=wheatprotein(:,1:6);
% alpha=0.01;
% u=lrt_env(X,Y,alpha)


function u=lrt_env(X,Y,alpha)

[n r]=size(Y);

stat0=env(X,Y,r);


for i=0:r-1
%     i
        stat=env(X,Y,i);
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