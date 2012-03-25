function WInit=get_Init4henv(F,X,Y,u,dataParameter)


n=dataParameter.n;
ng=dataParameter.ng;
p=dataParameter.p;
r=dataParameter.r;
sigY=dataParameter.sigY;
sigRes=dataParameter.sigRes;

tmp=zeros(r,r);
for i=1:p
    tmp=tmp+ng(i)/n*sigRes(:,:,i);
end

[V1 D]=eig(tmp);
[V2 D]=eig(sigY);
V=[V1 V2];

crit=nchoosek(2*r,u);

if crit<=50

    Ys=zeros(crit,r,u);
    [Ys(1:crit,:,:) ind]=get_combeig(V,2*r,u);
    imax = size(Ys,1);
    m = size(Ys,2);
    nn = size(Ys,3);
    
    minIndex = 1;
    minValue = F(reshape(Ys(minIndex,:,:),m,nn));

    for i=2:imax
        if (F(reshape(Ys(i,:,:),m,nn)) < minValue)
            minIndex = i;
            minValue = F(reshape(Ys(minIndex,:,:),m,nn));
        end
    end

    % the Y given the minimum value is our guess
    WInit=reshape(Ys(minIndex,:,:),m,nn);


else
    
    initset=zeros(1,u+1);
    initset(1:u)=1:u;
    initset(u+1)=1;
    iniValue = F(V(:,initset(1:u)));
    
    for rep=1:5
        for i=1:u+3
            for j=1:2*r
                if sum(j ==initset(2:u))==0
                    initset(1)=j;
                    temp=F(V(:,initset(1:u)));
                    if (temp < iniValue)
                        initset(u+1) = j;
                        iniValue = temp;
                    end
                end
            end %for j=1:2*r
            initset(1:u) = initset(2:(u+1));
            initset(u+1)=initset(1);
        end % end for i=1:u
    end % end for rep=1:3

    WInit=V(:,initset(1:u));
    
end

