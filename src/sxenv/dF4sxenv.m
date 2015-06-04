%% dF4sxenv
% The first derivative of the objective function for computing the envelope
% subspace for the scaled predictor envelope model.

%% Syntax
%         df = dF4sxenv(R, DataParameter)
% 
%% Input
%
% *R*: A p by u semi orthogonal matrix, 0 < u <= p.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *df*: An p by u matrix containing the value of the derivative function
% evaluated at R.

%% Description
%
% The objective function is derived in Section 2.2 of Cook and Su (2015) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF4sxenv(R, DataParameter)

n = DataParameter.n;
sigXcY = DataParameter.sigXcY;
invSigX = DataParameter.invSigX;
Lambda = DataParameter.Lambda;

a = 2 * Lambda \ sigXcY / Lambda * R / (R' / Lambda * sigXcY / Lambda * R);

b = 2 * Lambda * invSigX * Lambda * R / (R' * Lambda * invSigX * Lambda * R);

df = n * (a + b);