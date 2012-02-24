function Wguess=startv(X,Y,u)

%----- Get initial value for envelope model, using the algorithm in partial envelope paper --------
%----- Here we only use the MLE of beta since we assume that the dimension
%for outer envelope is r.
%----- if all the combination of r choose u is small(<=20), we find out 
%the combination that minimizes F, if that number is large, then we do it 
%iteratively: we pick up the first u eigenvalues of Sigma_est, fix all of 
%them except the first one and then search for the first one that minimizes 
%F, then we record that vector and fix all of them except the second, 
%search for the second and record the vector, that goes on and on until the 
%last one. We do it for 5 rounds and use the final set as our starting 
%value. 


global sigres;
global FParameters;
% sigma = FParameters.sigma;
% sigmag = FParameters.sigmag;

r=size(Y,2);
[n p]=size(X);

XC=center(X);
YC=center(Y);
[betfm sigres]=fit_OLS(X,Y);
betfm=betfm';
[V D]=eig(sigres);

crit=nchoosek(r,u);

if crit<=50
%if u==1
%     setaux(Y,X,0,r,V,D);
%     Ys=zeros(crit+5,r,u);
    Ys=zeros(crit,r,u);
    [Ys(1:crit,:,:) ind]=get_combeig(V,r,u);
%     Ys(crit+1,:,:) = getSIR(u,r);
%     Ys(crit+2,:,:) = getSAVE(u,r);
%     Ys(crit+3,:,:) = getDR(u,r);
%     Ys((crit+3+1):end,:,:) = get_more(sigma,sigmag,r,0,u);
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
    Wguess=reshape(Ys(minIndex,:,:),m,nn);
    bini=Wguess*pinv(Wguess'*Wguess)*Wguess'*betfm';
    sigini=(YC-XC*bini')'*(YC-XC*bini')/n;
        

    m=1;
    while (m*p<u)
        m=m+1;
    end
    m=m+2;
    spini=zeros(r,m*p);
    for i=1:m
        spini(:,(i-1)*p+1:i*p)=sigini^(i-1)*bini;
    end
    rk=min(m*p,r);
    
    [Ul S V]=svd(spini);
    [Ss In]=sort(diag(S(1:rk,1:rk)),'descend');
    Wguess=Ul(:,1:u);

else
    
    initset=zeros(1,u+1);
    initset(1:u)=1:u;
    initset(u+1)=1;
    iniValue = F(V(:,initset(1:u)));
    
    for rep=1:5
        for i=1:u+3
            for j=1:r
                if sum(j ==initset(2:u))==0
                    initset(1)=j;
                    temp=F(V(:,initset(1:u)));
                    if (temp < iniValue)
                        initset(u+1) = j;
                        iniValue = temp;
                    end
                end
            end %for j=1:r
            initset(1:u) = initset(2:(u+1));
            initset(u+1)=initset(1);
        end % end for i=1:u
    end % end for rep=1:3

    Wguess=V(:,initset(1:u));
    bini=Wguess*pinv(Wguess'*Wguess)*Wguess'*betfm';
    sigini=(YC-XC*bini')'*(YC-XC*bini')/n;
        

    m=1;
    while (m*p<u)
        m=m+1;
    end
    m=m+2;
    spini=zeros(r,m*p);
    for i=1:m
        spini(:,(i-1)*p+1:i*p)=sigini^(i-1)*bini;
    end
    rk=min(m*p,r);
    
    [Ul S V]=svd(spini);
    [Ss In]=sort(diag(S(1:rk,1:rk)),'descend');
    Wguess=Ul(:,1:u);
    
end

