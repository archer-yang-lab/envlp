%% Expan
% Compute the expansion matrix of dimension r.

%% Syntax
%         E = Expan(r)
%
%% Input
%
% *r*: Dimension of the expansion matrix.  A positive integer.
% 
%% Output
%
% *E*: Expansion matrix of dimension r.  E is an r ^ 2 by r(r + 1) / 2 matrix.
 
%% Description
% 
% The contraction and expansion matrices are links between the "vec" 
% operator and "vech" operator: for an r by r symmetric matrix A, 
% vech(A) = Contr(r) * vec(A), and vec(A) = Expan(r) * vech(A).  
% The "vec" operator stacks the matrix A into an r ^ 2 by 1 vector
% columnwise.  The "vech" operator stacks the lower triangle or the upper
% triangle of a symmetric matrix into an r(r + 1) / 2 vector. For more details
% of "vec", "vech", contraction and expansion matrix, refer to Henderson
% and Searle (1979).

function E = Expan(r)

% For future reference, the (i, j)th (j >= i) element of a matrix A, corresponds to
% the r * (i - 1) + j th element in vec(A), but corresponds to the
% (2 * r - i) / 2 * (i - 1) + j th element in vech(A). If j < i, it corresponds to the
% (2 * r - j) / 2 * (j - 1) + i th element in vech(A), but r * (i - 1) + j th element in
% vec(A).

    E = zeros(r ^ 2,r * (r + 1) / 2);
    
    for i = 1 : r
        for j = 1 : r
            if (j >= i)
                E(r * (i - 1) + j, (2 * r - i) / 2 * (i - 1) + j) = 1;
            else
                E(r * (i - 1) + j, (2 * r - j) / 2 * (j - 1) + i) = 1;
            end
        end
    end
    
    