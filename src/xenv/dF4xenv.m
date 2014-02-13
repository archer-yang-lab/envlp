%% dF4xenv
% The first derivative of the objective function for computing the envelope
% subspace for the reduction on X.

%% Syntax
%         df = dF4xenv(R, DataParameter)
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
% *df*: A p by u matrix containing the value of the derivative function
% evaluated at R.

%% Description
%
% The objective function is derived in Section 4.5.1 of Cook et al. (2013) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF4xenv(R, DataParameter)

n = DataParameter.n;
sigXcY = DataParameter.sigXcY;
invSigX = DataParameter.invSigX;

a = 2 * sigXcY * R / (R' * sigXcY * R);

b = 2 * invSigX * R / (R' * invSigX * R);

df = n * (a + b);