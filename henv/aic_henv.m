%% aic_henv
% Select the dimension of the envelope subspace using Akaike information
% criterion for the heteroscedastic envelope model.

%% Usage
% u=aic_henv(X,Y)
%
% Input
%
% * X: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% Output
%
% * u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace for the heteroscedastic envelope model. 

function u=aic_henv(X,Y)

[n r]=size(Y);
    
stat=henv(X,Y,r);
ic=-2*stat.l+2*stat.np;
u=r;


for i=0:r-1

        stat=henv(X,Y,i);
        temp=-2*stat.l+2*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end
