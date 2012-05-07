%% bic_senv
% Select the dimension of the scaled envelope subspace using Bayesian
% information criterion.

%% Syntax
% u=bic_senv(X,Y)
% u=bic_senv(X,Y,Opts)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors and n 
% is the number of observations. The predictors can be univariate or 
% multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.
% * Opts: The optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used.
%
% Output
%
% * u: Dimension of the inner envelope. An integer between 0 and r.
% 
%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the scaled envelope subspace. 

%% Example
%
% load('T9-12.txt')
% Y=T9_12(:,4:7);
% X=T9_12(:,1:3);
% u=bic_senv(X,Y)

function u=bic_senv(X,Y,Opts)

[n r]=size(Y);
    
ModelOutput=senv(X,Y,r,Opts);
ic=-2*ModelOutput.l+log(n)*ModelOutput.np;
u=r;


for i=0:r-1

        ModelOutput=senv(X,Y,i,Opts);
        temp=-2*ModelOutput.l+log(n)*ModelOutput.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
        
end