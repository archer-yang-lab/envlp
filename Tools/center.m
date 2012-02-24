function XC = center(X);
% XC = center(X)
% Centers the matrix X
%==========================
n = size(X,1);

XC=X-ones(n,1)*mean(X);

