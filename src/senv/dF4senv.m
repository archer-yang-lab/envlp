%% dF4senv
% First derivative of the objective function for computing the envelope
% subspace in the scaled envelope model.

%% Syntax
%         df = dF4senv(R, DataParameter)
% 
%% Input
%
% *R*: An r by u semi-orthogonal matrix, 0 < u <= r.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *df*: The first derivative of the objective function for computing the
%  envelope subspace.  An r by u matrix.

%% Description
%
% This first derivative of F4senv obtained by matrix calculus calculations.

function df = dF4senv(R, DataParameter)

n = DataParameter.n;
sigRes = DataParameter.sigRes;
Lambda = DataParameter.Lambda;
invsigY = DataParameter.invsigY;
invLambda = diag(1./diag(Lambda));

a = 2 * invLambda * sigRes * invLambda * R / (R' * invLambda * sigRes * invLambda * R);

b = 2 * Lambda * invsigY * Lambda * R / (R' * Lambda * invsigY * Lambda * R);

df = n * (a + b);