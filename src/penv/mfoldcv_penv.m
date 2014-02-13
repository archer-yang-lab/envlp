%% mfoldcv_penv
% Select the dimension of the partial envelope subspace using m-fold cross validation.

%% Syntax
%         u = mfoldcv_penv(X, Y, m)
%         u = mfoldcv_penv(X, Y, m, Opts)

%% Input
%
% *X*: A list containing the value of X1 and X2.
% 
% * X.X1: Predictors of main interest. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X.X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.
%
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables.
%
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.

%% Output
% 
%  *u*: The dimension of the partial envelope subspace selected by m-fold cross
%  validation. An integer between 0 and r.

%% Description
% This function implements m-fold cross validation to select the dimension
% of the partial envelope space, based on prediction performance.  For each u, the
% data is partitioned into m parts, each part is in turn used for testing 
% for the prediction performance while the rest m-1 parts are used for 
% training.  The dimension is selected as the one that minimizes the average 
% prediction errors. As Y is multivariate, the identity inner product is 
% used for computing the prediction errors.

%% Example
% 
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X.X1 = fiberpaper(:, 7);
%         X.X2 = fiberpaper(:, 5 : 6);
%         u = mfoldcv_penv(X, Y, 5)

function u = mfoldcv_penv(X, Y, m, Opts)

if nargin < 3
    error('Inputs: X, Y, and m should be specified!');
elseif nargin == 3
    Opts = [];
end

X.X1 = double(X.X1);
X.X2 = double(X.X2);
Y = double(Y);

X1 = X.X1;
X2 = X.X2;

[n, r] = size(Y);


Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

tempInd = min(floor((m - 1) * n / m) - 1, r);
PreErr = zeros(m, tempInd + 1);

for j = 0 : tempInd
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(j + 1) '\n']);
    end
    
    for i = 1 : m

        index = true(n, 1);
        index((floor((i - 1) * n / m) + 1) : ceil(i * n / m)) = 0;
        tempX.X1 = X1(index, :);
        tempX.X2 = X2(index, :);
        tempY = Y(index, :);
        ModelTemp = penv(tempX, tempY, j);
        
        testX.X1 = X1(logical(1 - index), :);
        testX.X2 = X2(logical(1 - index), :);
        testY = Y(logical(1 - index), :);
        testN = size(testY, 1);
        resi = testY - ones(testN, 1) * ModelTemp.alpha' - testX.X1 * ModelTemp.beta1' - testX.X2 * ModelTemp.beta2';
        PreErr(i, j + 1) = sqrt(trace(resi * resi') / testN);
        
    end
end


[~, ind] = min(mean(PreErr));
u = ind - 1;

