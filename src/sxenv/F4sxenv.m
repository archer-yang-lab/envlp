%% F4sxenv
% Objective function for computing the envelope subspace for the scaled 
% predictor envelope model.

%% Syntax
%         f = F4sxenv(R, DataParameter)
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
% *f*: A scalar containing the value of the objective function evaluated at R.

%% Description
%
% The objective function is derived in Section 2.2 of Cook and Su (2015)
%  using maximum likelihood estimation. The columns of the semi-orthogonal 
% matrix that minimizes this function span the estimated envelope subspace.

function f = F4sxenv(R, DataParameter)

n = DataParameter.n;
sigXcY = DataParameter.sigXcY;
invSigX = DataParameter.invSigX;
Lambda = DataParameter.Lambda;


eigtem = eig(R' / Lambda * sigXcY / Lambda * R);
a = sum(log(eigtem(eigtem > 0)));

eigtem0 = eig(R' * Lambda * invSigX * Lambda * R);
b = sum(log(eigtem0(eigtem0 > 0)));

f = n * (a + b);