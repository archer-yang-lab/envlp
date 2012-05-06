%% Lmatrix
% Extract the 2nd to the last diagonal element of a matrix into a vector.

%% Syntax
% L=Lmatrix(r)
%
% Input
%
% * r: The dimension of the matrix being extracted.  The matrix should be
% an r by r matrix.
%
% Output
%
% * L: An r-1 dimensional vector that contains all the diagonal elements
% but the first one of the matrix.

%% Description
% Let A be an r by r matrix, and vec be the vector operator, then 
% Lmatrix(r)*vec(A) will give the 2nd to the rth diagonal elements of A,
% arranged in a column vector. 

function L=Lmatrix(r)

L=zeros(r^2,r-1);
for i=1:r-1
    L(i*r+i+1,i)=1;
end