%% F4ienv
% Objective function for computing the inner envelope subspace.

%% Syntax
%         f = F4ienv(R, DataParameter)
% 
%% Input
%
% *R*: An r by u semi orthogonal matrix, 0 < u <= p.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *f*: A scalar containing the value of the objective function evaluated at R.

%% Description
%
% The objective function is derived in Section 3.3 in Su and Cook (2012) by
%  using maximum likelihood estimation. The columns of the semi-orthogonal 
% matrix that minimizes this function span the estimated inner envelope subspace.


function f = F4ienv(R, DataParameter)

u = size(R, 2);

r = DataParameter.r;
p = DataParameter.p;
n = DataParameter.n;
sigRes = DataParameter.sigRes;
sigFit = DataParameter.sigFit;
logDetSigRes = DataParameter.logDetSigRes;
invsigRes = DataParameter.invsigRes;

eigtem = eig(R' * sigRes * R);
a = log(prod(eigtem(eigtem > 0)));

eigtem0 = eig(R' * invsigRes * R);
b = log(prod(eigtem0(eigtem0 > 0)));

R0 = grams(nulbasis(R'));
[~, D] = eig(eye(r - u) / (R0' * sigRes * R0) * R0' * sigFit * R0);
lambdas = sort(diag(D), 'descend');
logl = log(lambdas + 1);
c = sum(logl((p - u + 1) : end));

f = n * r * (1 + log(2 * pi)) + n * (logDetSigRes + a + b + c);
