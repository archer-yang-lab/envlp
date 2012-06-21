%% get_Init4envmean
% Starting value for the envelope subspace in estimating the multiavriate mean.

%% Syntax
%         WInit = get_Init4envmean(F, u, DataParameter)

%% Input
%
% *F*: Objective function to get the envelope subspace.
% 
% *u*: Dimension of the envelope. An integer between 1 and p - 1.
% 
% *DataParameter*: A list containing commonly used statistics computed from
% the data.
%
%% Output
%
% *WInit*: The initial estimate of the orthogonal basis of the envelope
% subspace. A p by u orthogonal matrix.

%% Description
% We compute the eigenvectors for the estimated error covariance matrix,
% and get p vectors.  Then we get all the combinations  
% of u vectors out of the p vectors. If the number of p choose u is 
% small(<=50), we search over all the combinations and find out the one 
% that minimizes the objective function F. If that number is large, then we
% do it iteratively: we pick up any u eigenvectors, fix all of them except 
% the first one. Then we search over all the vectors orthogonal to the 
% fixed ones, and record the one that minimizes F. Next, we fix the first 
% u eigenvectors again but this time search for the second one, then we record 
% the vector. This goes on and on until the last one. We do it for 3 rounds 
% and use the final set as our starting value. 

%% Reference
% The codes are implemented based on the algorithm in Section 3.5 of Su and 
% Cook (2011).

function WInit = get_Init4envmean(F, u, DataParameter)

n = DataParameter.n;
p = DataParameter.p;
sigX = DataParameter.sigX;


[V D] = eig(sigX);

crit = nchoosek(p, u);

if crit <= 50
%if u==1
%     setaux(Y,X,0,r,V,D);
%     Ys=zeros(crit+5,r,u);
    Ys = zeros(crit, p, u);
    [Ys(1 : crit, :, :) ind] = get_combeig(V, p, u);
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
    WInit = reshape(Ys(minIndex, :, :), m, nn);

else
    
    initset = zeros(1, u + 1);
    initset(1 : u) = 1 : u;
    initset(u + 1) = 1;
    iniValue = F(V(:, initset(1 : u)));
    
    for rep = 1 : 5
        for i = 1 : u + 3
            for j = 1 : p
                if sum(j == initset(2 : u)) == 0
                    initset(1) = j;
                    temp = F(V(:, initset(1 : u)));
                    if (temp < iniValue)
                        initset(u + 1) = j;
                        iniValue = temp;
                    end
                end
            end 
            initset(1 : u) = initset(2 : (u + 1));
            initset(u + 1) = initset(1);
        end 
    end 

    WInit = V(:, initset(1 : u));
    
end

