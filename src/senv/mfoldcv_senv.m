%% mfoldcv_senv
% Select the dimension of the scaled envelope subspace 
% using m-fold cross validation.

%% Syntax
%         u = mfoldcv_senv(X, Y, m)
%         u = mfoldcv_senv(X, Y, m, Opts)

%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations. 
% 
% *Y*: Responses. An n by r matrix, r is the number of responses. The 
% number of the responses should be greater than the number of the predictors.
% And they must be continuous variables.  
%
% *m*: A positive integer that is used to indicate m-fold cross validation.
% 
% *Opts*: A list containing the optional input parameters. If one or
% several (even all) fields are not defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.
% * Opts.rep: Number of replicates for scales. This option imposes special 
% structure on scaling parameters. For example, if Opts.rep = [3 4], this 
% means that the first three responses have the same scale and the next 
% four responses share a different scale. The elements of this vector should 
% sum to r. If not specified, the default is [], then all responses will be
% scaled differently. If all responses have the same scale, input [r], then 
% the regular envelope will be applied to the data.
% The input should be a row vector.

%% Output
% 
%  *u*: The dimension of the scaled envelope subspace selected by m-fold cross
%  validation. An integer between 0 and r.

%% Description
% This function implements m-fold cross validation to select the dimension
% of the scaled envelope space, based on prediction performance.  For each u, the
% data is partitioned into m parts, each part is in turn used for testing 
% for the prediction performance while the rest m-1 parts are used for 
% training.  The dimension is selected as the one that minimizes the average 
% prediction errors. As Y is multivariate, the identity inner product is 
% used for computing the prediction errors.

%% Example
% 
%         load('sales.txt')
%         Y = sales(:, 4 : 7);
%         X = sales(:, 1 : 3);
%         u = mfoldcv_senv(X, Y, 5)

function u = mfoldcv_senv(X, Y, m, Opts)

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

        index = true(n, 1);
        index((floor((i - 1) * n / m) + 1) : ceil(i * n / m)) = 0;
        tempX = X(index, :);
        tempY = Y(index, :);
        ModelTemp = senv(tempX, tempY, j);
        
        testX = X(logical(1 - index), :);
        testY = Y(logical(1 - index), :);
        testN = size(testX, 1);
        resi = testY - ones(testN, 1) * ModelTemp.alpha' - testX * ModelTemp.beta';
        PreErr(i, j + 1) = sqrt(trace(resi * resi') / testN);
        
    end
end


[~, ind] = min(mean(PreErr));
u = ind - 1;

