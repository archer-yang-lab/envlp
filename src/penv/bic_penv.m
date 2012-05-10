%% bic_penv
% Select the dimension of the partial envelope subspace using Bayesian information
% criterion.

%% Syntax
%         u = bic_penv(X, Y)
%         u = bic_penv(X, Y, Opts)
%
%% Input
%
% *X*: A list containing the value of X1 and X2.
% 
% * X.X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
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
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the partial envelope subspace. 

%% Example
%         load T7-7.dat
%         Y = T7_7(:, 1 : 4);
%         Xtemp = T7_7(:, 5 : 7);
%         X.X1 = Xtemp(:, 3);
%         X.X2 = Xtemp(:, 1 : 2);
%         u = bic_penv(X, Y)

function u = bic_penv(X, Y, Opts)

if nargin < 2
    error('Inputs: X and Y should be specified!');
elseif nargin == 2
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n r] = size(Y);
    
ModelOutput = penv(X, Y, r, Opts);
ic = -2 * ModelOutput.l + log(n) * ModelOutput.np;
u = r;


for i = 0 : r - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = penv(X, Y, i, Opts);
    temp = - 2 * ModelOutput.l + log(n) * ModelOutput.np;
    
    if temp < ic
        u = i;
        ic = temp;
    end
    
end
