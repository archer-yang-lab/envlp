%% dF4env
% The first derivative of the objective funtion for computing the envelope
% subspace.

%% Syntax
% df = dF4env(R, dataParameter)
% 
% Input
%
% * R: An r by u semi orthogonal matrix, 0<u<=r.
% * dataParameter: A structure that contains the statistics calculated from
% the data.
%
% Output
%
% * df: An r by u matrix containing the value of the derivative function
% evaluated at R.

%% Description
%
% The objective function is derived in Section 4.3 in Cook et al. (2010) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF4env(R, dataParameter)

sigRes = dataParameter.sigRes;
sigY = dataParameter.sigY;

a = 2 * sigRes * R * inv(R' * sigRes * R);

temp = inv(sigY);

b = 2 * temp * R * inv(R' * temp * R);

df = a + b;