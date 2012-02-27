%% Contr
% Compute the contraction matrix of dimension r.

%% Usage
% C = Contr(r)
%
% Input
%
% * r: Dimension of the contraction matrix.  A positive integer.
% 
% Output
%
% * C: Contraction matrix of dimension r.  C is an r(r+1)/2 by r^2 matrix.
 
%% Description
% 
% The contraction and expansion matrices are links between the "vec" 
% operator and "vech"operator: for an r by r symmetric matrix A, 
% vech(A)=Contr(r)vec(A), and vec(A)=Expan(r)vech(A).  
% The "vec" operator stacks the matrix A into an r^2 by 1 vector
% columnwise.  The "vech" operator stacks the lower triangle or the upper
% triangle of a symmetric matrix into an r(r+1)/2 vector. For more details
% of "vec", "vech", contraction and expansion matrix, refer to Henderson
% and Searle (1979).

function C = Contr(r)

% For future reference, the (i,j)th (j>=i) element of a matrix A, corresponds to
% the r*(i-1)+j th element in vec(A), but corresponds to the
% (2*r-i)/2*(i-1)+j th element in vech(A). If j<i, it corresponds to the
% (2*r-j)/2*(j-1)+i th element in vech(A), but r*(i-1)+j th element in
% vec(A).

    C=zeros(r*(r+1)/2,r^2);
    
    for i=1:r
        for j=1:r
            if (j==i)
                C((2*r-i)/2*(i-1)+j,r*(i-1)+j)=1;
            elseif (j>i)
                C((2*r-i)/2*(i-1)+j,r*(i-1)+j)=1/2;
            else
                C((2*r-j)/2*(j-1)+i,r*(i-1)+j)=1/2;
            end
        end
    end