%% lrt_henv
% Select the dimension of the envelope subspace using likelihood ratio
% testing for the heteroscedastic envelope model.

%% Syntax
%         u = lrt_henv(X, Y, alpha)
%         u = lrt_henv(X, Y, alpha, Opts)
%
%% Input
%
% *X*: Group indicators. A matrix with n rows.  X can only have p unique
%  rows, where p is the number of groups. For example, if there 
% are two groups, X can only have 2 different kinds of rows, such as (0, 1)
% and (1, 0), or (1, 0, 10) and (0, 5, 6).  The number of columns is not
% restricted, as long as X only has p unique rows.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be greater than p.
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
% *u*: Dimension of the envelope. An integer between 0 and r.

%% Description
% This function implements the likelihood ratio testing procedure to select
% the dimension of the envelope subspace in heteroscedastic envelope model,
% with pre-specified significance level $$\alpha$. 

%% Example
% 
%         load waterstrider.mat
%         u = lrt_henv(X, Y, 0.01)

function u = lrt_henv(X, Y, alpha, Opts)

if (nargin < 3)
    error('Inputs: X, Y and alpha should be specified!');
elseif (nargin == 3)
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n r] = size(Y);

ModelOutput0 = henv(X, Y, r, Opts);


for i = 0 : r - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = henv(X, Y, i, Opts);
    chisq = - 2 * (ModelOutput.l - ModelOutput0.l);
    df = ModelOutput0.np - ModelOutput.np;
    
    if chi2cdf(chisq, df) < (1 - alpha)
        u = i;
        break;
    end
    
end

if (i == r - 1) && chi2cdf(chisq, df) > (1 - alpha)
    u = r;
end