%% Kpd
% Compute the communication matrix Kpd.

%% Syntax
%         k = Kpd(p, d)
%
%% Input
%
% *p*, *d*: two positive integers represent the dimension parameters
% for the communication matrix.
%
%% Output
%
% *k*: The communication matrix Kpd. An p * d by p * d matrix.

%% Description
% For a p by d matrix A, vec(A') = Kpd * vec(A), and Kpd is called a
% communication matrix.

%% Reference
% The codes are implemented based on Definition 3.1 in Magnus and Neudecker
% (1979).

function k = Kpd(p, d)

k = zeros(p * d, p * d);

for i = 1 : p
    for j = 1 : d
        Hij = zeros(p, d);
        Hij(i, j) = 1;
        k = k + kron(Hij, Hij');
    end
end