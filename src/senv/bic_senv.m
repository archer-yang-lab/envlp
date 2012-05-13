%% bic_senv
% Select the dimension of the scaled envelope subspace using Bayesian
% information criterion.

%% Syntax
%         u = bic_senv(X, Y)
%         u = bic_senv(X, Y, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n 
% is the number of observations. The predictors can be univariate or 
% multivariate, discrete or continuous.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.
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
% *u*: Dimension of the inner envelope. An integer between 0 and r.
% 
%% Description
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the scaled envelope subspace. 

%% Example
%
%         load('sales.txt')
%         Y = sales(:, 4 : 7);
%         X = sales(:, 1 : 3);
%         u = bic_senv(X, Y)

function u = bic_senv(X, Y, Opts)

if nargin < 2
    error('Inputs: X and Y should be specified!');
elseif nargin == 2
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n r] = size(Y);
    
ModelOutput = senv(X, Y, r, Opts);
ic = - 2 * ModelOutput.l + log(n) * ModelOutput.np;
u = r;


for i = 0 : r - 1
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
    ModelOutput = senv(X, Y, i, Opts);
    temp = - 2 * ModelOutput.l + log(n) * ModelOutput.np;
    
    if temp < ic
        u = i;
        ic = temp;
    end
    
end
