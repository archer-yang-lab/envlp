%% mfoldcv_henv
% Use m-fold cross validation to select the dimension of the envelope subspace 
% for heteroscedastic envelope model.

%% Syntax
%         u = mfoldcv_henv(X, Y, m)
%         u = mfoldcv_henv(X, Y, m, Opts)

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
% continuous variables.
%
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag to print out dimension selection process, 
% logical 0 or 1. Default value: 0.
% * Opts.table: Flag to tabulate the results, which contains cross 
% validation error for each u.  Logical 0 or 1. Default value: 0.

%% Output
% 
%  *u*: The dimension of the envelope subspace selected by m-fold cross
%  validation.  An integer between 0 and r.

%% Description
% This function implements m-fold cross validation to select the dimension
% of the envelope space, based on prediction performance.  For each u, the
% data is partitioned into m parts, each part is in turn used for testing 
% for the prediction performance while the rest m-1 parts are used for 
% training.  The dimension is selected as the one that minimizes the average 
% prediction errors. As Y is multivariate, the identity inner product is 
% used for computing the prediction errors.

%% Example
% 
%         load waterstrider.mat
%         u = mfoldcv_henv(X, Y, 5)

function u = mfoldcv_henv(X, Y, m, Opts)

if nargin < 3
    error('Inputs: X, Y, and m should be specified!');
elseif nargin == 3
    Opts = [];
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

X = double(X);
Y = double(Y);

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n, p] = size(X);

tempInd = min(floor((m - 1) * n / m) - 1, p);
PreErr = zeros(1, tempInd + 1);

for j = 0 : tempInd
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(j + 1) '\n']);
    end
    
    for i = 1 : m

        index = true(n, 1);
        index((floor((i - 1) * n / m) + 1) : floor(i * n / m)) = 0;
        tempX = X(index, :);
        tempY = Y(index, :);
        ModelTemp = henv(tempX, tempY, j);
        
        testX = X(logical(1 - index), :);
        testY = Y(logical(1 - index), :);
        testN = size(testX, 1);
        sqe = 0;
        for k = 1 : testN
            pred = predict_henv(ModelTemp, testX(k, :)', 'estimation');
            sqe = sqe + (testY(k, :) - pred.value') * (testY(k, :)' - pred.value);
        end
        PreErr(j + 1) = PreErr(j + 1) + sqe;
        
    end
    
    PreErr(j + 1) = sqrt(PreErr(j + 1) / n);
end


[~, ind] = min(PreErr);
u = ind - 1;

if tableFlag == 1
    
    fprintf('\n u      CV error      \n');
    fprintf('------------------------\n');
    for i = 0 : tempInd
        fprintf('%2d %12.3f\n', i, PreErr(i + 1));
    end
    fprintf('------------------------\n');
    
end