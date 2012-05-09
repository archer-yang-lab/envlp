%% fit_OLS
% Multivariate linear regression. 

%% Syntax
% ModelOutput = fit_OLS(X, Y)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
%
% Output
%
% ModelOutput: A list that contains the maximum likelihood estimators of
% regression coefficients and error covariance matrix. 
% 
% * ModelOutput.betaOLS: An r by p matrix containing estimate of the regression coefficients $$\beta$.
% * ModelOutput.SigmaOLS: An r by r matrix containing estimate of the error covariance matrix.
% * ModelOutput.alpha: An r by 1 vector containing estimate of the intercept.
% * ModelOutput.n: The number of observations in the data.  A positive integer.

%% Description
% In a multivariate linear model, Y and X follows the following
% relationship: $$Y=\alpha+\beta X+\varepsilon$, where $$\varepsilon$
% contains the errors.  This function performs the ordinary least squares
% fit to the inputs, and returns the estimates of $$\beta$ and the
% covariance matrix of $$\varepsilon$.

%% Example
%
% load wheatprotein.txt
% X = wheatprotein(:, 8);
% Y = wheatprotein(:, 1 : 6);
% ModelOutput = fit_OLS(X, Y)

function ModelOutput = fit_OLS(X, Y)

n = length(X);
XC = center(X);
YC = center(Y);
mY = mean(Y)';
mX = mean(X)';
PX = XC * inv(XC' * XC) * XC';
ModelOutput.betaOLS = YC' * XC * inv(XC' * XC);
ModelOutput.SigmaOLS = YC' * (eye(n) - PX) * YC / n;
ModelOutput.alpha = mY - ModelOutput.betaOLS * mX;
ModelOutput.n = n;