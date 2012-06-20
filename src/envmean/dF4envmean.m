%% dF4envmean
% The first derivative of the objective function for computing the envelope
% subspace.

%% Syntax
%         df = dF4envmean(R, DataParameter)
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
% The objective function is derived by maximum likelihood estimation. This
% function is the derivative of the objective function. 

function df = dF4envmean(R, DataParameter)

n = DataParameter.n;
sX = DataParameter.sX;
sigX = DataParameter.sigX;

a = 2 * sigX * R * inv(R' * sigX * R);

temp = inv(sX);

b = 2 * temp * R * inv(R' * temp * R);

df = n * (a + b);