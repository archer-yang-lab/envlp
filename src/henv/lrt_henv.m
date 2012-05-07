%% lrt_henv
% Select the dimension of the envelope subspace using likelihood ratio
% testing for the heteroscedastic envelope model.

%% Syntax
% u=lrt_henv(X,Y,alpha)
% u=lrt_henv(X,Y,alpha,Opts)
%
% Input
%
% * X: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% * alpha: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
% * Opts: A list containing the optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used. 
%
% Output
%
% * u: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the envelope subspace in heteroscedastic envelope model,
% with prespecified significance level $$\alpha$. 

%% Example
% 
% load waterstrider.mat
% u=lrt_henv(X,Y,0.01)

function u=lrt_henv(X,Y,alpha,Opts)

if (nargin < 3)
    error('Inputs: X, Y and alpha should be specified!');
elseif (nargin==3)
    Opts=[];
end

[n r]=size(Y);

ModelOutput0=henv(X,Y,r,Opts);


for i=0:r-1

        ModelOutput=henv(X,Y,i,Opts);
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