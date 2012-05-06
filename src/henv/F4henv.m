%% F4henv
% Objective funtion for computing the envelope subspace in heteroscedastic
% envelope model.

%% Usage
% f = F4henv(R,dataParameter)
% 
% Input
%
% * R: An r by u semi orthogonal matrix, 0<u<=r.
% * dataParameter: A structure that contains the statistics calculated from
% the data.
%
% Output
%
% * f: A scalar containing the value of the objective function evaluated at
% R.

%% Description
%
% The objective function is derived in Section 2.2 of Su and Cook (2012)
%  using maximum likelihood estimation. The columns of the semi-orthogonal 
% matrix that minimizes this function span the estimated envelope subspace
% in the heteroscedastic envelope model.

function f = F4henv(R,dataParameter)

p=dataParameter.p;
n=dataParameter.n;
ng=dataParameter.ng;
sigRes=dataParameter.sigRes;
sigY=dataParameter.sigY;

f=0;
for i=1:p
    eigtem = eig(R' * sigRes(:,:,i) * R);
    a = log(prod(eigtem(eigtem>0)));
    f = f+ng(i)/n*a;  
end


eigtem0=eig(R'*inv(sigY)*R);
b=log(prod(eigtem0(eigtem0>0)));

f=f+b;