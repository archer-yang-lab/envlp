function dataParameter=make_parameter4henv(X,Y)

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

