%% aic_henv
% Select the dimension of the envelope subspace using Akaike information
% criterion for the heteroscedastic envelope model.

%% Syntax
% u=aic_henv(X,Y)
% u=aic_henv(X,Y,opts)
%
% Input
%
% * X: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
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
% the dimension of the envelope subspace for the heteroscedastic envelope model. 

%% Example
% 
% load waterstrider.mat
% u=aic_henv(X,Y)

function u=aic_henv(X,Y,opts)

if (nargin < 2)
    error('Inputs: X, Y should be specified!');
elseif (nargin==2)
    opts=[];
end

[n r]=size(Y);
    
ModelOutput=henv(X,Y,r,opts);
ic=-2*ModelOutput.l+2*ModelOutput.np;
u=r;


for i=0:r-1

        ModelOutput=henv(X,Y,i,opts);
        temp=-2*ModelOutput.l+2*ModelOutput.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
        
end
