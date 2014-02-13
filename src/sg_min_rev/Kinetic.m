function K = Kinetic(n)

% Generates a square tridiagonal matrix, with the main diagonal elements
% all being 2, subdiagonal and superdiagonal elements all being -1. 
	K = 2*eye(n) - diag(ones(1,n-1),-1) - diag(ones(1,n-1),1);
