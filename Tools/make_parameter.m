function dataParameter=make_parameter(X,Y,method)

if (strcmp(method,'env'))
    
    
    [n p]=size(X);
    r=size(Y,2);

    XC=center(X);
    YC=center(Y);
    [betaOLS sigRes]=fit_OLS(X,Y);

    dataParameter.n=n;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.XC=XC;
    dataParameter.YC=YC;
    dataParameter.mX=mean(X)';
    dataParameter.mY=mean(Y)';
    dataParameter.sigX=cov(X,1);
    dataParameter.sigY=cov(Y,1);
    dataParameter.sigRes=sigRes;
    dataParameter.betaOLS=betaOLS;
   
    
elseif (strcmp(method,'senv'))
    
    
    [n p]=size(X);
    r=size(Y,2);

    XC=center(X);
    YC=center(Y);
    [betaOLS sigRes]=fit_OLS(X,Y);

    dataParameter.n=n;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.mX=mean(X)';
    dataParameter.mY=mean(Y)';
    dataParameter.sigX=cov(X,1);
    dataParameter.sigY=cov(Y,1);
    dataParameter.sigRes=sigRes;
    dataParameter.betaOLS=betaOLS;
    
    
elseif (strcmp(method,'ienv'))
    
    
    [n p]=size(X);
    r=size(Y,2);

    XC=center(X);
    YC=center(Y);
    [betaOLS sigRes]=fit_OLS(X,Y);
    sigY=cov(Y,1);
    sigFit=sigY-sigRes;

    dataParameter.n=n;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.XC=XC;
    dataParameter.YC=YC;
    dataParameter.mX=mean(X)';
    dataParameter.mY=mean(Y)';
    dataParameter.sigX=cov(X,1);
    dataParameter.sigY=sigY;
    dataParameter.sigRes=sigRes;
    dataParameter.sigFit=sigFit;
    dataParameter.betaOLS=betaOLS;
    
    
elseif (strcmp(method,'henv'))
    
    
    [n r]=size(Y);

    [Xs ind]=sortrows(X);
    p=1;
    temp=Xs(1,:);
    ng=0;
    for i=1:n
        if prod(double(Xs(i,:)==temp))==0            
            temp=Xs(i,:);
            ng(p)=i-1-sum(ng(1:p-1));
            ncum(p)=i-1;
            p=p+1;
        end
    end           
    ncum(p)=n;
    ng(p)=n-sum(ng(1:p-1));

    sigRes=zeros(r,r,p);
    mYg=zeros(r,p);

    Ys=Y(ind,:);

    for i=1:p
        if i>1
            sigRes(:,:,i)=cov(Ys(ncum(i-1)+1:ncum(i),:),1);
            mYg(:,i)=mean(Ys(ncum(i-1)+1:ncum(i),:))';
        else
            sigRes(:,:,i)=cov(Ys(1:ncum(i),:),1);
            mYg(:,i)=mean(Ys(1:ncum(i),:))';
        end
    end



    dataParameter.n=n;
    dataParameter.ng=ng;
    dataParameter.ncum=ncum;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.ind=ind;
    dataParameter.mY=mean(Y)';
    dataParameter.mYg=mYg;
    dataParameter.sigY=cov(Y,1);
    dataParameter.sigRes=sigRes;
    
    

end


