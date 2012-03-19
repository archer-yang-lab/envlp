%% F4ienv
% Objective funtion for computing the inner envelope subspace

%% Usage
% f = F4ienv(R,dataParameter)
% 
% Input
%
% * R: An r by u semi orthogonal matrix, 0<u<=p.
% * dataParameter: A structure that contains the statistics calculated form
% the data.
%
% Output
%
% * f: A scalar containing the value of the objective function evaluated at R.

%% Description
%
% The objective function is derived in Section 3.3 in Su and Cook (2012) by
%  using maximum likelihood estimation. The columns of the semi-orthogonal 
% matrix that minimizes this function span the estimated inner envelope subspace.


function f = F4ienv(R,dataParameter)

u=size(R,2);

sigRes=dataParameter.sigRes;
sigY=dataParameter.sigY;
sigFit=sigY-sigRes;
p=dataParameter.p;



    eigtem=eig(R'*sigRes*R);
    a=log(prod(eigtem(eigtem>0)));

    eigtem0=eig(R'*inv(sigRes)*R);
    b=log(prod(eigtem0(eigtem0>0)));

    R0=grams(nulbasis(R'));
    [V D]=eig(inv(R0'*sigRes*R0)*R0'*sigFit*R0);
    lambdas=sort(diag(D),'descend');
    logl=log(lambdas+1);
    c=sum(logl((p-u+1):end));    

    f=a+b+c;
