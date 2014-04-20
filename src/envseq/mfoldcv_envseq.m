%% mfoldcv_envseq
% Select the dimension of the envelope subspace using m-fold cross
% validation for envelope model using a sequential
% algorithm.

%% Syntax
%         u = mfoldcv_envseq(X, Y, m)
%         u = mfoldcv_envseq(X, Y, m, Opts)

%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations. 
% 
% *Y*: Responses. An n by r matrix, r is the number of
% responses. The responses must be continuous variables.  
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
%  validation.

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
%         load Rohwer  
%         X = Rohwer(:, 4 : 5); 
%         Y = Rohwer(:, 1 : 3);
%         m = 5;
%         u = mfoldcv_envseq(X, Y, m)

function u = mfoldcv_envseq(X, Y, m, Opts)

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

[n, r] = size(Y);

tempInd = min(floor((m - 1) * n / m) - 1, r);
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
        ModelTemp = envseq(tempX, tempY, j);
        
        testX = X(logical(1 - index), :);
        testY = Y(logical(1 - index), :);
        testN = size(testX, 1);
        resi = testY - ones(testN, 1) * ModelTemp.alpha' - testX * ModelTemp.beta';
        PreErr(j + 1) = PreErr(j + 1) + trace(resi * resi');
        
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