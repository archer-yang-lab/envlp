%% objfun_spls
% Objective function for computing a-star in the scaled SIMPLS algorithm.

%% Syntax
%         f = objfun_spls(a, Gamma1, Gamma2, DataParameter)
% 
%% Input
%
% *a*: A scalar between 0 and 1 containing the initial value of a.
% 
% *Gamma1*: A p by u matrix containing an orthogonal basis of the envelope
% subspace in the last iteration.
% 
% *Gamma2*: A p by u matrix containing an orthogonal basis of the envelope
% subspace in the current iteration.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *f*: A scalar containing the value of the objective function evaluated at a.

%% Description
%
% The objective function is introduced in Section 3 of Cook and Su (2015)
% under the context of scaled SIMPLS algorithm (SPLS algorithm). 

function f = objfun_spls(a, Gamma1, Gamma2, DataParameter)

    n = DataParameter.n;
    invSigX = DataParameter.invSigX;
    sigXcY = DataParameter.sigXcY;
    Gamma = grams(a * Gamma1 + (1 - a) * Gamma2);
    d = DataParameter.d;
    rep = DataParameter.rep;

    C = arrayfun(@(x, y) repmat(x, [1 y]), [1 d], rep, 'UniformOutput', false);
    Ld = cell2mat(C);
    Lambda = diag(Ld);
    eigtem = eig(Gamma' * Lambda * invSigX * Lambda * Gamma);
    a = sum(log(eigtem(eigtem > 0))); 
    
    eigtem2 = eig(Gamma' / Lambda * sigXcY / Lambda * Gamma);
    b = sum(log(eigtem2(eigtem2 > 0)));
    
    f = n * a / 2 + n * b / 2;