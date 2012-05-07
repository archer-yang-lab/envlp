%% bic_penv
% Select the dimension of the partial envelope subspace using Bayesian information
% criterion.

%% Syntax
% u=bic_penv(X1,X2,Y)
% u=bic_penv(X1,X2,Y,Opts)
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
% * Opts: A list containing the optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used. 
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

function u=bic_penv(X1,X2,Y,Opts)

if (nargin < 3)
    error('Inputs: X1, X2 and Y should be specified!');
elseif (nargin==3)
    Opts=[];
end

[n r]=size(Y);
    
ModelOutput=penv(X1,X2,Y,r,Opts);
ic=-2*ModelOutput.l+log(n)*ModelOutput.np;
u=r;


for i=0:r-1

        ModelOutput=penv(X1,X2,Y,i,Opts);
        temp=-2*ModelOutput.l+log(n)*ModelOutput.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
        
end
