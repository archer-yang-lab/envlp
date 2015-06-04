%% objfun
% Objective function for computing the scales in the scaled envelope model.

%% Syntax
%         f = objfun(d, Gamma, DataParameter)
% 
%% Input
%
% *d*: An r - 1 dimensional column vector containing the scales for the 2nd
% to the rth responses.  All the entries in d are positive.
% 
% *Gamma*: A r by u semi-orthogonal matrix that spans the envelope subspace
% or the estimated envelope subspace.
% 
% *DataParameter*: A structure that contains the statistics calculated form
% the data.
%
%% Output
%
% *f*: A scalar containing the value of the objective function evaluated at d.

%% Description
%
% The objective function is derived in Section 4.1 of Cook and Su (2013)
%  using maximum likelihood estimation. 

function f = objfun(d, Gamma, DataParameter)

n = DataParameter.n;
sigRes = DataParameter.sigRes;
rep = DataParameter.rep;
invsigY = DataParameter.invsigY;

    C = arrayfun(@(x, y) repmat(x, [1 y]), [1 d], rep, 'UniformOutput', false);
    Ld = cell2mat(C);
    Lambda = diag(Ld);
    invLambda = diag(1./Ld);
    eigtem = eig(Gamma' * Lambda * invsigY * Lambda * Gamma);
    a = sum(log(eigtem(eigtem > 0))); 
    eigtem2 = eig(Gamma' * invLambda * sigRes * invLambda * Gamma);
    b = sum(log(eigtem2(eigtem2 > 0)));
    
    f = n * a / 2 + n * b / 2;