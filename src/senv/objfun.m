%% objfun
% Objective funtion for computing the scales in the scaled envelope model.

%% Syntax
% f = objfun(d,Gamma,DataParameter)
% 
% Input
%
% * d: An r-1 dimensional column vector containing the scales for the 2nd
% to the rth responses.  All the entries in d are positive.
% * Gamma: A r by u semi-orthogomal matrix that spans the envelope subspace
% or the estimated envelope subspace.
% * DataParameter: A structure that contains the statistics calculated form
% the data.
%
% Output
%
% * f: A scalar containing the value of the objective function evaluated at d.

%% Description
%
% The objective function is derived in Section 4.1 of Su and Cook (2012)
%  using maximum likelihood estimation. 

function f=objfun(d,Gamma,DataParameter)

n=DataParameter.n;
sigY=DataParameter.sigY;
sigRes=DataParameter.sigRes;



    Lambda=diag([1 d]);
    eigtem = eig(Gamma'*Lambda*inv(sigY)*Lambda*Gamma);
    a= log(prod(eigtem(eigtem>0))); 
    eigtem2 = eig(Gamma'*inv(Lambda)*sigRes*inv(Lambda)*Gamma);
    b= log(prod(eigtem2(eigtem2>0)));
    
    f = n*a/2+n*b/2;