%% bic_ienv
% Select the dimension of the inner envelope subspace using Bayesian
% information criterion.

%% Syntax
% u = bic_ienv(X, Y)
% u = bic_ienv(X, Y, Opts)
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
% * u: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.

%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the inner envelope subspace. 

%% Example
%
% load irisf.mat
% u = bic_ienv(X, Y)

function u = bic_ienv(X, Y, Opts)

if nargin < 2
    error('Inputs: X, Y should be specified!');
elseif nargin == 2
    Opts = [];
end

[n r] = size(Y);
p = size(X, 2);
    
ModelOutput = env(X, Y, r, Opts);
ic = - 2 * ModelOutput.l + log(n) * ModelOutput.np;
u = r;


for i = 0 : p

        ModelOutput = ienv(X, Y, i, Opts);
        temp = - 2 * ModelOutput.l + log(n) * ModelOutput.np;
        
        if temp < ic
           u = i;
           ic = temp;
        end
end
