%% lrt_xenv
% 
% Use likelihood ratio testing to select the dimension of the envelope
% subspace for the reduction on X.

%% Syntax
%         u = lrt_xenv(X, Y, alpha)
%         u = lrt_xenv(X, Y, alpha, Opts)
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
% *Opts*: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the envelope. An integer between 0 and p.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the envelope subspace for the reduction on X, with
% pre-specified significance level $$\alpha$. 

%% Example
%
%         load wheatprotein.txt
%         X = wheatprotein(:, 1 : 6);
%         Y = wheatprotein(:, 7);
%         alpha = 0.01;
%         u = lrt_xenv(X, Y, alpha)

function u = lrt_xenv(X, Y, alpha, Opts)

if nargin < 3
    error('Inputs: X, Y and alpha should be specified!');
elseif nargin == 3
    Opts = [];
end

if (alpha < 0 || alpha > 1)
    error('alpha should be between [0, 1]!');
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n p] = size(X);

ModelOutput0 = xenv(X, Y, p, Opts);


for i = 0 : p - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = xenv(X, Y, i, Opts);
    chisq = - 2 * (ModelOutput.l - ModelOutput0.l);
    df = ModelOutput0.np - ModelOutput.np;
    
    if chi2cdf(chisq, df) < (1 - alpha)
        u = i;
        break;
    end
    
end

if i == p - 1 && chi2cdf(chisq, df) > (1 - alpha)
    u = p;
end