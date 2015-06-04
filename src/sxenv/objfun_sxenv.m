%% objfun_sxenv
% Objective function for computing the scales for the scaled predictor 
% envelope model.

%% Syntax
%         f = objfun_sxenv(d, Gamma, DataParameter)
% 
%% Input
%
% *d*: A (p - 1) by 1 vector containing the value of the scales (the second 
% to the last diagonal element in $\Lambda$.
% 
% *Gamma*: A p by u matrix containing an orthogonal basis of the envelope
% subspace.
% 
% *DataParameter*: A structure that contains the statistics calculated from
% the data.
%
%% Output
%
% *f*: A scalar containing the value of the objective function evaluated at d.

%% Description
%
% The objective function is derived in Section 2.2 of Cook and Su (2015) by
%  using maximum likelihood estimation. This function is the derivative of
%  the objective function.

function f = objfun_sxenv(d, Gamma, DataParameter)

    n = DataParameter.n;
    invSigX = DataParameter.invSigX;
    sigXcY = DataParameter.sigXcY;
    rep = DataParameter.rep;

    C = arrayfun(@(x, y) repmat(x, [1 y]), [1 d], rep, 'UniformOutput', false);
    Ld = cell2mat(C);
    Lambda = diag(Ld);
    eigtem = eig(Gamma' * Lambda * invSigX * Lambda * Gamma);
    a = sum(log(eigtem(eigtem > 0))); 
    
    eigtem2 = eig(Gamma' / Lambda * sigXcY / Lambda * Gamma);
    b = sum(log(eigtem2(eigtem2 > 0)));
    
    f = n * a / 2 + n * b / 2;