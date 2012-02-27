%% dF
% The derivative of the objective funtion for computing the envelope subspace

%% Usage
% df = dF ( R )
% 
% Input
%
% * R: An r by u semi orthogonal matrix, 0<u<=r.
%
% Output
%
% * df: An r by u containing the value of the derivative function evaluated at R.

%% Description
%
% The objective function is derived in Section 4.3 in Cook et al. (2010) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function df = dF(R)

global sigY;
global sigres;

a=sigres*R*inv(R'*sigres*R);

temp=inv(sigY);

b=temp*R*inv(R'*temp*R);

df=a+b;