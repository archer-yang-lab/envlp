%% lrt_ienv
% Select the dimension of the inner envelope subspace using likelihood
% ratio testing.

%% Syntax
% u = lrt_ienv(X, Y, alpha)
% u = lrt_ienv(X, Y, alpha, Opts)
%
%% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% * alpha: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
% * Opts: The optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used.
%
%% Output
%
% * u: Dimension of the inner envelope. An integer between 0 and p or equal
% to r.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the inner envelope subspace, with prespecified significance 
% level $$\alpha$.

%% Example
%
% load irisf.mat
% 
% alpha = 0.01;
% u = lrt_ienv(X, Y, alpha)


function u = lrt_ienv(X, Y, alpha, Opts)

if nargin < 3
    error('Inputs: X, Y and alpha should be specified!');
elseif nargin == 3
    Opts = [];
end

[n r] = size(Y);
p = size(X, 2);

ModelOutput0 = env(X, Y, r, Opts);


for i = 1 : (p + 1)

        ModelOutput = ienv(X, Y, p + 1 - i, Opts);
        chisq = - 2 * (ModelOutput.l - ModelOutput0.l);
        df = ModelOutput0.np - ModelOutput.np;
        
        if chi2cdf(chisq, df) < (1 - alpha)
            u = p + 1 - i;
            break;
        end
        
end

if i == p + 1 && chi2cdf(chisq, df) > (1 - alpha)
    u = r;
    warning('No inner envelope model is selected, fit with the standard multivariate linear model.');
end