%% startv
% Starting value for the envelope subspace.

%% Usage
% WInit=get_Init(X,Y,u,dataParameter)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. 
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations.  
% * u: Dimension of the envelope. An integer between 1 and r-1.
%
% Output
%
% * WInit: The initial estimate of the orthogonal basis of the envelope
% subspace. An r by u orthogonal matrix.

%% Description
% We compute the eigenvectors for the covariance matrices of Y and the 
% estimated errors, and get 2r vectors.  Then we get all the combinations 
% of u vectors out of the 2r vectors. If the number of 2r choose u is 
% small(<=50), we search over all the combinations and find out the one 
% that minimizes the objective function F. If that number is large, then we
% do it iteratively: we pick up any u eigenvectors, fix all of them except 
% the first one. Then we search over all the vectors orthogonal to the 
% fixed ones, and record the one that minimizes F. Next, we fix the first 
% u eigenvectors again but this time search for second one, then we record 
% the vector. This goes on and on until the last one. We do it for 5 rounds 
% and use the final set as our starting value. 

%% Reference
% The codes is implemented based on the algorithm in Section 3.5 of Su and 
% Cook (2011).

function WInit=get_Init(F,X,Y,u,dataParameter)


% global sigres;
% global FParameters;
% sigma = FParameters.sigma;
% sigmag = FParameters.sigmag;

n=dataParameter.n;
p=dataParameter.p;
r=dataParameter.r;
XC=dataParameter.XC;
YC=dataParameter.YC;
sigY=dataParameter.sigY;
sigRes=dataParameter.sigRes;
betaOLS=dataParameter.betaOLS;


[V1 D]=eig(sigRes);
[V2 D]=eig(sigY);
V=[V1 V2];

crit=nchoosek(2*r,u);

if crit<=50
%if u==1
%     setaux(Y,X,0,r,V,D);
%     Ys=zeros(crit+5,r,u);
    Ys=zeros(crit,r,u);
    [Ys(1:crit,:,:) ind]=get_combeig(V,2*r,u);
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
    bini=Wguess*pinv(Wguess'*Wguess)*Wguess'*betaOLS;
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
    WInit=Ul(:,1:u);

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
    bini=Wguess*pinv(Wguess'*Wguess)*Wguess'*betaOLS;
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
    WInit=Ul(:,1:u);
    
end

