%% dF4xenv
% The first derivative of the objective funtion for computing the envelope
% subspace for the reduction on X.

%% Syntax
% df = dF4xenv(R, DataParameter)
% 
% Input
%
% * R: An r by u semi orthogonal matrix, 0 < u <= p.
% * DataParameter: A structure that contains the statistics calculated from
% the data.
%
% Output
%
% * df: An p by u matrix containing the value of the derivative function
% evaluated at R.

%% Description
%
% The objective function is derived in Section 4.5.1 of Cook et al. (2012) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF4xenv(R, DataParameter)

sigXcY = DataParameter.sigXcY;
sigX = DataParameter.sigX;
invSigX = DataParameter.invSigX;

a = 2 * sigXcY * R * inv(R' * sigXcY * R);

b = 2 * invSigX * R * inv(R' * invSigX * R);

df = a + b;