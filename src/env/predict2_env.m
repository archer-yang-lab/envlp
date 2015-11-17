%% predict2_env
% Perform estimation or prediction under the envelope model through partial envelope model.

%% Syntax
%         PredictOutput = predict2_env(X, Y, u, Xnew, infType)
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
% *u*: The dimension of the constructed partial envelope model.  An integer
% between from 0 to r.
% 
% *Xnew*: The value of X with which to estimate or predict Y.  A p by 1
% vector.
% 
% *infType*: A string of characters indicating the inference type,
% the choices can be 'estimation' or 'prediction'.
% 
%% Output
%
% *PredictOutput*: A list containing the results of the inference.
%
% * PredictOutput.value: The fitted value or the prediction value evaluated at
% Xnew. An r by 1 vector.
% * PredictOutput.covMatrix: The covariance matrix of PredictOutput.value. An r by r
% matrix.
% * PredictOutput.SE: The standard error of elements in PredictOutput.value.  An r
% by 1 vector. 

%% Description
% 
% This function evaluates the envelope model at new value Xnew.  It can
% perform estimation: find the fitted value when X = Xnew, or prediction:
% predict Y when X = Xnew.  The covariance matrix and the standard errors are
% also provided.  Compared to predict_env, this function performs 
% prediction through partial envelope model, which can be more accurate if
% the partial envelope is of smaller dimension and contains less variant
% material information.

%% References
% 
% # The codes are implemented based on the following reference: R.D. Cook 
% (2013) ``Lecture Notes on Envelope Models and Methods.'' School of
% Statistics, University of Minnesota, Minneapolis. 

%% Example
%
%         load fiberpaper.dat
%         Y = fiberpaper(:, 1 : 4);
%         X = fiberpaper(:, [5 6 7]);
%         alpha = 0.01;
%         u = lrt_env(X, Y, alpha);
%         ModelOutput = env(X, Y, u);
%         Xnew = X(10, :)';
%         p1 = predict_env(ModelOutput, Xnew, 'estimation')
%         p1.value
%         p1.SE
%         u = lrt_predict2_env(X, Y, 0.01, Xnew)
%         p2 = predict2_env(X, Y, 1, Xnew, 'estimation')
%         p2.value
%         p2.SE
%         p1.SE./p2.SE

function PredictOutput = predict2_env(X, Y, u, Xnew, infType)

if nargin < 5
    error('Inputs: X, Y, u, Xnew and infType should be specified!');
end

if ~strcmp(infType, 'estimation') && ~strcmp(infType, 'prediction')
    error('Inference type can only be estimation or prediction.');
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

X0 = grams(nulbasis(Xnew'));
A = grams([Xnew X0]);
Ainv = inv(A);
Z = X * (Ainv)';

Xtemp.X1 = Z(:, 1);
Xtemp.X2 = Z(:, 2 : end);

mtemp = penv(Xtemp, Y, u);
Ntemp.X1 = Ainv(1, :) * Xnew;
Ntemp.X2 = Ainv(2 : end, :) * Xnew;
PredictOutput = predict_penv(mtemp, Ntemp, infType);

