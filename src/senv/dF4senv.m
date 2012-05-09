%% dF4senv
% First derivative of the objective funtion for computing the envelope
% subspace in the scaled envelope model.

%% Syntax
% f = dF4senv(R, DataParameter)
% 
% Input
%
% * R: An r by u semi-orthogonal matrix, 0 < u <= p.
% * DataParameter: A structure that contains the statistics calculated from
% the data.
%
% Output
%
% * dF: The first derivative of the objective function for computing the
%  envelope subspace.  An r by u matrix.

%% Description
%
% This first derivative of F4senv obtained by matrix calculus calculations.

function df = dF4senv(R, DataParameter)

sigRes = DataParameter.sigRes;
sigY = DataParameter.sigY;
Lambda = DataParameter.Lambda;

a = 2 * inv(Lambda) * sigRes * inv(Lambda) * R * inv(R' * inv(Lambda) * sigRes * inv(Lambda) * R);

temp = inv(sigY);

b = 2 * Lambda * temp * Lambda * R * inv(R' * Lambda * temp * Lambda * R);

df = a + b;