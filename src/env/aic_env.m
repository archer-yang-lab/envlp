%% aic_env
% Select the dimension of the envelope subspace using Akaike information
% criterion.

%% Syntax
% u=aic_env(X,Y)
% u=aic_env(X,Y,opts)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% * opts: A list containing the optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used. 
%
% Output
%
% * u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace. 

%% Example
% load wheatprotein.txt
% X=wheatprotein(:,8);
% Y=wheatprotein(:,1:6);
% u=aic_env(X,Y)

function u=aic_env(X,Y,opts)

if (nargin < 2)
    error('Inputs: X, Y should be specified!');
elseif (nargin==2)
    opts=[];
end

[n r]=size(Y);
    
stat=env(X,Y,r,opts);
ic=-2*stat.l+2*stat.np;
u=r;


for i=0:r-1
    
        stat=env(X,Y,i,opts);
        temp=-2*stat.l+2*stat.np;
        
        if (temp<ic)
           u=i;
           ic=temp;
        end
        
end
