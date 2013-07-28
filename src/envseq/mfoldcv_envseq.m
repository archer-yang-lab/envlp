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
% responses. The 
% number of the responses should be greater than the number of the predictors.
% And they must be continuous variables.  
%
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out dimension seletion process, 
% logical 0 or 1. Default value: 0.

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
% prediction errors. If Y is multivariate, the identity inner product is 
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

X = double(X);
Y = double(Y);

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n, r]=size(Y);

tempInd = min(floor((m - 1) * n / m) - 1, r);
PreErr = zeros(m, tempInd + 1);

for j = 0 : tempInd
    
    if printFlag == 1
        fprintf(['Current dimension ' int2str(j + 1) '\n']);
    end
    
    for i = 1 : m

        index = logical([ones(n, 1)]);
        index((floor((i - 1) * n / m) + 1) : ceil(i * n / m)) = 0;
        tempX = X(index, :);
        tempY = Y(index, :);
        ModelTemp = envseq(tempX, tempY, j);
        
        testX = X(logical(1 - index), :);
        testY = Y(logical(1 - index), :);
        testN = size(testX, 1);
        resi = testY - ones(testN, 1) * ModelTemp.alpha' - testX * ModelTemp.beta';
        PreErr(i, j + 1) = sqrt(trace(resi * resi') / testN);
        
    end
end


[minErr ind] = min(mean(PreErr));
u = ind - 1;

