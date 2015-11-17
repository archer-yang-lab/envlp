%% bic_predict2_env
% Select the dimension of the constructed partial envelope subspace using Bayesian information
% criterion.
%% Syntax
%         u = bic_predict2_env(X, Y, Xnew)
%         u = bic_predict2_env(X, Y, Xnew, Opts)
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
% *Xnew*: The value of X with which to estimate or predict Y.  A p by 1
% vector.
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
% * Opts.table: Flag to tabulate the results, which contains AIC and log
% likelihood for each u. Logical 0 or 1. Default value: 0.
% 
%% Output
%
% *u*: Dimension of the constructed partial envelope. An integer between 0 and r.

%% Description
% 
% This function implements the Bayesian information criteria (BIC) to select
% the dimension of the partial envelope model, which is constructed for prediction based on envelope model. 

%% References
% 
% # The codes are implemented based on the following reference: R.D. Cook 
% (2013) ``Lecture Notes on Envelope Models and Methods.'' School of
% Statistics, University of Minnesota, Minneapolis. 

%% Example
%
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X = fiberpaper(:, [7 5 6]);
%         Xnew = X(10, :)';
%         Opts.table = 1;
%         u = bic_predict2_env(X, Y, Xnew, Opts)


function u = bic_predict2_env(X, Y, Xnew, Opts)

if nargin < 3
    error('Inputs: X, Y and Xnew should be specified!');
elseif nargin == 3
    Opts = [];
end

[n, p] = size(X);
n1 = size(Y, 1);
if n ~= n1
    error('The number of observations in X and Y should be equal!');
end

[s1, s2] = size(Xnew);

if s1 ~= p || s2 ~= 1
    error('Xnew must be a p by 1 vector');
end

if p == 1
    error('This method does not apply to p = 1.')
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


X0 = grams(nulbasis(Xnew'));
A = grams([Xnew X0]);
Ainv = inv(A);
Z = X * (Ainv)';

Xtemp.X1 = Z(:, 1);
Xtemp.X2 = Z(:, 2 : end);

u = bic_penv(Xtemp, Y, Opts);


