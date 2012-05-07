%% lrt_env
% Select the dimension of the envelope subspace using likelihood ratio
% testing.

%% Syntax
% u = lrt_env(X, Y, alpha)
% u = lrt_env(X, Y, alpha, opts)
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
% * opts: The optional input parameter. If one or several (even all) 
% fields are not defined, the default settings (see make_opts documentation) 
% are used.
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
% X = wheatprotein(:, 8);
% Y = wheatprotein(:, 1:6);
% alpha = 0.01;
% u = lrt_env(X, Y, alpha)


function u = lrt_env(X, Y, alpha, opts)

if nargin < 3
    error('Inputs: X, Y and alpha should be specified!');
elseif nargin == 3
    opts = [];
end

[n r] = size(Y);

stat0 = env(X, Y, r, opts);

for i = 0 : r-1
	stat = env(X, Y, i, opts);
	chisq = - 2 * (stat.l - stat0.l);
	df = stat0.np - stat.np;
	if (chi2cdf(chisq, df) < (1 - alpha))
	    u = i;
	    break;
	end
end

if i == r-1 && chi2cdf(chisq, df) > (1 - alpha)
    u = r;
end