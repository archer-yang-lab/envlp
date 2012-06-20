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

