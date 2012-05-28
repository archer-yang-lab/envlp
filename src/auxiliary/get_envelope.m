%% get_envelope
% Construct the envelope subspace using a sequential algorithm.

%% Syntax
%         W = get_envelope(S, M, u)
%
%% Input
%
% *S*: An r by p matrix whose columns span the subspace, the rank of S
% cannot be greater than u. 
% 
% *M*: An r by r positive semi-definite matrix.
%
% *u*: Dimension of the envelope. An integer between 0 and r.
%

%% Output
%
% *W*: An r by u semi-orthogonal matrix that spans the M-envelope of
% span(S).

%% Description
% This function constructs the M-envelope of span(S) using a sequential
% algorithm similar to partial least squares.

%% Reference
% # The codes are implemented based on the algorithm in the lecture notes
% of Cook (2012).


function W = get_envelope(S, M, u)

if nargin < 3
    error('Inputs: S, M and u should be specified!');
end

r = size(M, 1);
r1 = size(S, 1);
p = rank(S);

if r ~= r1
    error('The size of S is not valid.');
end

if p > u
    error('The rank of S cannot be greater than u.');
end

if sum(eig(M) < 0) > 0
    error('M must be semi-positive definite.');
end

U = S * S';

[w, D] = eig(U);
[DM ind] = max(diag(D));
W = w(:, ind);

if u > 1
    
    for i = 2 : u
        
        Eu = grams(M * W);
        QEu = eye(r) - Eu * inv(Eu' * Eu) * Eu';
        [w, D] = eigs(QEu * U * QEu, 1);
        W = [W w];
        
    end
    
end


