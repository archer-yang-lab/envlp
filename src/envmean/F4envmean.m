%% F4envmean
% Objective function for computing the envelope subspace.

%% Syntax
%         f = F4envmean(R, DataParameter)
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
% *f*: A scalar containing the value of the objective function evaluated at
% R. 

%% Description
%
% The objective function is derived by maximum likelihood estimation.
% The columns of the semi-orthogonal matrix that minimizes this function
% span the estimated envelope subspace.  

function f = F4envmean(R, DataParameter)

n = DataParameter.n;
p = DataParameter.p;
sigY = DataParameter.sigY;
logDetSY = DataParameter.logDetSY;
invsY = DataParameter.invsY;

eigtem = eig(R' * sigY * R);
a = sum(log(eigtem(eigtem > 0)));

eigtem0 = eig(R' * invsY * R);
b = sum(log(eigtem0(eigtem0 > 0)));

f = n * p * (1 + log(2 * pi)) + n * (a + b + logDetSY);