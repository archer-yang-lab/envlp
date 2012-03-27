%% fit_OLS
% Multivariate linear regression. 

%% Usage
% [betaOLS SigmaOLS]=fit_OLS(X,Y)
%
% Input
%
% * X: Predictors, an n by p matrix, p is the number of predictors.  The predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
%
% Output
%
% * betaOLS: An r by p matrix containing estimate of the regression coefficients $$\beta$.
% * SigmaOLS: An r by r matrix containing estimate of the error covariance matrix.

%% Description
% In a multivariate linear model, Y and X follows the following
% relationship: $$Y=\alpha+\beta X+\varepsilon$, where $$\varepsilon$
% contains the errors.  This function performs the ordinary least squares
% fit to the inputs, and returns the estimates of $$\beta$ and the
% covariance matrix of $$\varepsilon$.


function [betaOLS SigmaOLS]=fit_OLS(X,Y)

n=length(X);
XC=center(X);
YC=center(Y);
PX=XC*inv(XC'*XC)*XC';
betaOLS=YC'*XC*inv(XC'*XC);
SigmaOLS=YC'*(eye(n)-PX)*YC/n;