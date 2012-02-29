%% F
% Objective funtion for computing the envelope subspace

%% Usage
% f = F( R )
% 
% Input
%
% * R: An r by u semi orthogonal matrix, 0<u<=r.
%
% Output
%
% * f: A scalar containing the value of the objective function evaluated at R.

%% Description
%
% The objective function is derived in Section 4.3 in Cook et al. (2010) by
%  using maximum likelihood estimation. The columns of the semi-orthogonal 
% matrix that minimizes this function span the estimated envelope subspace.

function f = F(R)

global sigY;
global sigres;

eigtem=eig(R'*sigres*R);
a=log(prod(eigtem(eigtem>0)));

eigtem0=eig(R'*inv(sigY)*R);
b=log(prod(eigtem0(eigtem0>0)));

f=a+b;