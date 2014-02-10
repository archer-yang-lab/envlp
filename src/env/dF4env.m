%% dF4env
% The first derivative of the objective function for computing the envelope
% subspace.

%% Syntax
%         df = dF4env(R, DataParameter)
% 
%% Input
%
% *R*: An r by u semi orthogonal matrix, 0 < u <= r.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *df*: An r by u matrix containing the value of the derivative function
% evaluated at R.

%% Description
%
% The objective function is derived in Section 4.3 in Cook et al. (2010) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF4env(R, DataParameter)

n = DataParameter.n;
sigRes = DataParameter.sigRes;
invsigY = DataParameter.invsigY;

a = 2 * sigRes * R / (R' * sigRes * R);

b = 2 * invsigY * R / (R' * invsigY * R);

df = n * (a + b);