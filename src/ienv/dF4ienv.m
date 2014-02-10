%% dF4ienv
% First derivative of the objective function for computing the inner envelope
% subspace.

%% Syntax
%         df = dF4ienv(R, DataParameter)
% 
%% Input
%
% *R*: An r by u semi-orthogonal matrix, 0 < u <= p.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *df*: The first derivative of the objective function for computing the
% inner envelope subspace.  An r by u matrix.

%% Description
%
% This first derivative of F4ienv obtained by matrix calculus calculations.

function df = dF4ienv(R, DataParameter)

[r, u] = size(R);

sigRes = DataParameter.sigRes;
sigFit = DataParameter.sigFit;
p = DataParameter.p;
n = DataParameter.n;
invsigRes = DataParameter.invsigRes;


a = 2 * sigRes * R / (R' * sigRes * R) + 2 * invsigRes * R / (R' * invsigRes * R);

R0 = grams(nulbasis(R'));
temp1 = eye(r - u) / (R0' * sigRes * R0);
temp2 = R0' * sigFit * R0;
dzdg0 = kron(eye(r - u), temp1 * R0' * sigFit) ...
    + Kpd(r - u, r - u) * kron(temp1, R0' * sigFit) ...
    - kron(temp2 * temp1, temp1 * R0' * sigRes) ...
    - Kpd(r - u, r - u) * kron(temp1, temp2 * temp1 * R0' * sigRes);
dg0dg1 = - kron(R0', R) * Kpd(r, u);

[V0, D0] = eig(temp1);
MultiPlier = V0 * diag(sqrt(diag(D0))) * V0';
[Vtmp, Dtmp] = eig(MultiPlier * temp2 * MultiPlier);
[Ds, ~] = sort(diag(Dtmp), 'descend');
b = zeros(1, (r - u) ^ 2);


for i = p - u + 1 : r - u
    V = eye(r - u) / MultiPlier * Vtmp(:, i);
    U = MultiPlier * Vtmp(:, i);
    b = b + 1 / (1 + Ds(i)) * reshape(V * U', 1, (r - u) ^ 2);
end

b = b * dzdg0 * dg0dg1;

df = n * (a + reshape(b, r, u));
