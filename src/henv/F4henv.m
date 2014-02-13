%% F4henv
% Objective function for computing the envelope subspace in heteroscedastic
% envelope model.

%% Syntax
%         f = F4henv(R, DataParameter)
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
% *f*: A scalar containing the value of the objective function evaluated at
% R.

%% Description
%
% The objective function is derived in Section 2.2 of Su and Cook (2013)
%  using maximum likelihood estimation. The columns of the semi-orthogonal 
% matrix that minimizes this function span the estimated envelope subspace
% in the heteroscedastic envelope model.

function f = F4henv(R, DataParameter)

p = DataParameter.p;
r = DataParameter.r;
n = DataParameter.n;
ng = DataParameter.ng;
sigRes = DataParameter.sigRes;
sigY = DataParameter.sigY;
logDetSigY = DataParameter.logDetSigY;

f = 0;
for i = 1 : p
    eigtem = eig(R' * sigRes(:, :, i) * R);
    a = log(prod(eigtem(eigtem > 0)));
    f = f + ng(i) / n * a;  
end

eigtem0 = eig(R' / sigY * R);
b = log(prod(eigtem0(eigtem0 > 0)));

f = n * r * (1 + log(2 * pi)) + n * (f + b + logDetSigY);