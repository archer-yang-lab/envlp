%% lrt_penv
% Select the dimension of the partial envelope subspace using likelihood 
% ratio testing.

%% Syntax
% u=lrt_penv(X1,X2,Y,alpha)
% u=lrt_penv(X1,X2,Y,alpha,opts)
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
% * opts: A list containing the optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used. 
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

function u=lrt_penv(X1,X2,Y,alpha,opts)

if (nargin < 4)
    error('Inputs: X1, X2, Y and alpha should be specified!');
elseif (nargin==4)
    opts=[];
end

[n r]=size(Y);

ModelOutput0=penv(X1,X2,Y,r,opts);


for i=0:r-1

        ModelOutput=penv(X1,X2,Y,i,opts);
        chisq = -2*(ModelOutput.l-ModelOutput0.l);
        df=ModelOutput0.np-ModelOutput.np;
        
        if (chi2cdf(chisq,df) < (1-alpha))
            u=i;
            break;
        end
        
end

if (i== r-1) && chi2cdf(chisq,df) > (1-alpha)
    u=r;
end