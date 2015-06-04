%% lrt_ienv
% Select the dimension of the inner envelope subspace using likelihood
% ratio testing.

%% Syntax
%         u = lrt_ienv(X, Y, alpha)
%         u = lrt_ienv(X, Y, alpha, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
% 
% *alpha*: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
% 
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag to print out dimension selection process. 
% Logical 0 or 1. Default value: 0.
% * Opts.table: Flag to tabulate the results, which contains log 
% likelihood, test statistic, degrees of freedom and p-value for each test. 
% Logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the inner envelope. An integer between 0 and p.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the inner envelope subspace, with pre-specified significance 
% level $$\alpha$.

%% Example
%
%         load irisf.mat
% 
%         alpha = 0.01;
%         u = lrt_ienv(X, Y, alpha)


function u = lrt_ienv(X, Y, alpha, Opts)

if nargin < 3
    error('Inputs: X, Y and alpha should be specified!');
elseif nargin == 3
    Opts = [];
end

if (alpha < 0 || alpha > 1)
    error('alpha should be between [0, 1]!');
end

if isfield(Opts, 'table')
    if (Opts.table ~= 1)
        tableFlag = 0;
    else
        tableFlag = 1;
    end
else
    tableFlag = 0;
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

p = size(X, 2);
llik = zeros(p + 1, 1);
tstat = zeros(p + 1, 1);
df = zeros(p + 1, 1);
pv = zeros(p + 1, 1);

ModelOutput0 = ienv(X, Y, 0, Opts);
llik(1) = ModelOutput0.l;

for i = 1 : p
    
    if printFlag == 1
        fprintf(['\n Current dimension ' int2str(p + 1 - i)]);
    end
    
    ModelOutput = ienv(X, Y, p + 1 - i, Opts);
	llik(p + 1 - i) = ModelOutput.l;
    tstat(p + 1 - i) = - 2 * (ModelOutput.l - ModelOutput0.l);
	df(p + 1 - i) = ModelOutput0.paramNum - ModelOutput.paramNum;
	pv(p + 1 - i) = 1 - chi2cdf(tstat(p + 1 - i), df(p + 1 - i));
    if pv(p + 1 - i) > alpha
	    u = p + 1 - i;
	    break;
    end
    
end

if i == p && pv(p + 1 - i) < alpha
    u = 0;
    fprintf('No inner envelope model is selected, fit with the standard multivariate linear model.');
end

if tableFlag == 1
    
    fprintf('\n u      log likelihood      test statistic     degrees of freedom    p-value\n');
    fprintf('----------------------------------------------------------------------------------\n');
    for i = 1 : p
        fprintf('%2d %15.3f %18.3f %18d %18.3f\n', p + 1 - i, llik(p + 1 - i), tstat(p + 1 - i), df(p + 1 - i), pv(p + 1 - i));
    end
    fprintf('%2d %15.3f\n', 0, llik(1));    
    fprintf('----------------------------------------------------------------------------------\n');
    
end