%% get_Init
% Starting value for the envelope subspace.

%% Syntax
%         WInit = get_Init(F, u, DataParameter)
%
%% Input
%
% *F*: Objective function of the envelope subspace.
% 
% *u*: Dimension of the envelope. An integer between 1 and r - 1.
% 
% *DataParameter*: A list containing commonly used statistics computed from
% the data.
%
%% Output
%
% *WInit*: The initial estimate of the orthogonal basis of the envelope
% subspace. An r by u orthogonal matrix.

%% Description
% We compute the eigenvectors for the covariance matrices of Y and the 
% estimated errors, and get 2r vectors.  Then we get all the combinations 
% of u vectors out of the 2r vectors. If the number of 2r choose u is 
% small(<= 50), we search over all the combinations and find out the one 
% that minimizes the objective function F. If that number is large, then we
% do it iteratively: we pick up any u eigenvectors, fix all of them except 
% the first one. Then we search over all the vectors orthogonal to the 
% fixed ones, and record the one that minimizes F. Next, we fix the first 
% u eigenvectors again but this time search for the second one, then we record 
% the vector. This goes on and on until the last one. We do it for 5 rounds 
% and use the final set as our starting value. 

%% Reference
% The codes are implemented based on the algorithm in Section 3.5 of Su and 
% Cook (2011).

function WInit = get_Init(F, u, DataParameter)



n = DataParameter.n;
p = DataParameter.p;
r = DataParameter.r;
XC = DataParameter.XC;
YC = DataParameter.YC;
sigY = DataParameter.sigY;
sigRes = DataParameter.sigRes;
betaOLS = DataParameter.betaOLS;
sigYX = YC' * XC / n;
initflag = DataParameter.initflag;

[V1, ~] = eig(sigRes);
[V2, ~] = eig(sigY);
V = [V1 V2];

if r <= 10
    crit = nchoosek(2 * r, u);
else
    crit = 2 * r;
end

if u == 1 || (r <= 10 && crit <= 50)
%if u==1
%     setaux(Y,X,0,r,V,D);
%     Ys=zeros(crit+5,r,u);
    Ys = zeros(crit, r, u);
    [Ys(1 : crit, :, :), ~] = get_combeig(V, 2 * r, u);
%     Ys(crit+1,:,:) = getSIR(u,r);
%     Ys(crit+2,:,:) = getSAVE(u,r);
%     Ys(crit+3,:,:) = getDR(u,r);
%     Ys((crit+3+1):end,:,:) = get_more(sigma,sigmag,r,0,u);
    imax = size(Ys, 1);
    m = size(Ys, 2);
    nn = size(Ys, 3);
    
    minIndex = 1;
    minValue = F(reshape(Ys(minIndex, :, :), m, nn));

    for i = 2 : imax
        if (F(reshape(Ys(i, :, :), m, nn)) < minValue)
            minIndex = i;
            minValue = F(reshape(Ys(minIndex, :, :), m, nn));
        end
    end

    % the Y given the minimum value is our guess
    Wguess1 = reshape(Ys(minIndex, :, :), m, nn);
    bini = Wguess1 * pinv(Wguess1' * Wguess1) * Wguess1' * betaOLS;
    sigini = (YC - XC * bini')' * (YC - XC * bini') / n;
        
    Wguess1 = get_envelope(bini, sigini, u);


    [V, D1] = eig(sigRes);
    tmp1 = betaOLS * sigYX';
    [~, ind] = sort(diag(V' * tmp1 * V), 'descend');
    Wguess2 = V(:, ind(1:u));


    [V2, D2] = eig(sigY);
    [~, ind] = sort(diag(V2' * tmp1 * V2), 'descend');
    init1 = V2(:, ind(1:u));
    
    if (F(init1) < F(Wguess2))
        Wguess2 = init1;
    end
    
    if (initflag == 1)
        tmp2 = inv(sqrt(D1));
        [~, ind] = sort(diag(tmp2 * V' * tmp1 * V * tmp2), 'descend');
        init1 = V(:, ind(1:u))
        if (F(init1) < F(Wguess2))
            Wguess2 = init1;
        end
        
        tmp2 = inv(sqrt(D2));
        [~, ind] = sort(diag(tmp2 * V2' * tmp1 * V2 * tmp2), 'descend');
        init1 = V2(:, ind(1:u)) 
        if (F(init1) < F(Wguess2))
            Wguess2 = init1;
        end
    end
    
    if (F(Wguess1) < F(Wguess2))
        WInit = Wguess1;
    else
        WInit = Wguess2;
    end

else
    
    [V, D1] = eig(sigRes);
    tmp1 = betaOLS * sigYX';
    [~, ind] = sort(diag(V' * tmp1 * V), 'descend');
    Wguess2 = V(:, ind(1:u));


    [V2, D2] = eig(sigY);
    [~, ind] = sort(diag(V2' * tmp1 * V2), 'descend');
    init1 = V2(:, ind(1:u));
    
    if (F(init1) < F(Wguess2))
        Wguess2 = init1;
    end
    
    if (initflag == 1)
        tmp2 = inv(sqrt(D1));
        [~, ind] = sort(diag(tmp2 * V' * tmp1 * V * tmp2), 'descend');
        init1 = V(:, ind(1:u)); 
        if (F(init1) < F(Wguess2))
            Wguess2 = init1;
        end
        
        tmp2 = inv(sqrt(D2));
        [~, ind] = sort(diag(tmp2 * V2' * tmp1 * V2 * tmp2), 'descend');
        init1 = V2(:, ind(1:u)); 
        if (F(init1) < F(Wguess2))
            Wguess2 = init1;
        end
    end
    
    WInit = Wguess2;
    
end

