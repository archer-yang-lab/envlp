% -----Test envelope model-----
% Generate data
n=100;
r=10;
u=2;
p=2;

X=rand(n,p);
GG=rand(r,r);
GG=grams(GG);
G=GG(:,1:u);
G0=GG(:,u+1:end);
sigma=1;
sigma0=5;
ita=rand(u,p);
bet=G*ita*10;
Sigma=G*G'*sigma^2+G0*G0'*sigma0^2;
epsil=mvnrnd(zeros(1,r),Sigma,n);
Y=X*bet'+epsil;

% Fit and check
stat=env(X,Y,u);
eig(stat.Omega)
eig(stat.Omega0)
norm(stat.beta-bet)
[beta_OLS sigres]=fit_OLS(X,Y);
norm(beta_OLS-bet)
subspace(stat.Gamma,G)

alpha=0.01;
u=lrt_env(Y,X,alpha)
u=bic_env(Y,X)
u=aic_env(Y,X)



bootse=bstrp_OLS(X,Y,B)
bootse=bstrp_env(X,Y,B,u)

load wheatprotein.txt
X=wheatprotein(:,8);
Y=wheatprotein(:,1:6);
alpha=0.01;
u=lrt_env(Y,X,alpha)
stat=env(X,Y,u);
stat.Omega
eig(stat.Omega0)
stat.ratio

u=bic_env(Y,X)
u=aic_env(Y,X)


%-----Test partial envelope model-----

load T7-7.dat
Y=T7_7(:,1:4);
X=T7_7(:,5:7);
X1=X(:,3);
X2=X(:,1:2);
alpha=0.01;
u=lrt_penv(X1,X2,Y,alpha)
u=lrt_env(Y,X,0.01);
stat=env(X,Y,u);
stat=penv(X1,X2,Y,1)

bootse=bstrp_penv(X1,X2,Y,B,u)

X1=X(:,3);
X2=X(:,1:2);
u=lrt_penv(X1,X2,Y,0.01)
stat=penv(X1,X2,Y,u)

X1=X(:,2);
X2=X(:,[1 3]);
u=lrt_penv(X1,X2,Y,0.01)
stat=penv(X1,X2,Y,u)

X1=X(:,1);
X2=X(:,2:3);
u=lrt_penv(X1,X2,Y,0.01)
stat=penv(X1,X2,Y,u)


%-----Test inner envelope model-----


rand('state',11);
randn('state',11);
r=10;
p=8;
d=2;
n=200;
C1=grams(rand(r,r));
Gamma1=C1(:,1:d);
Gamma0=C1(:,d+1:end);
D=diag([1 5 10 50 100 500 1000 5000 10000 50000]);
Sigma=Gamma1*D(1:d,1:d)*Gamma1'+Gamma0*D(d+1:end,d+1:end)*Gamma0';
Mean=zeros(r,1);
eta1=[eye(d) zeros(d,p-d)]';
eta2=[zeros(p-d,d) randn(p-d,p-d)]';

B=grams(rand(r-d,p-d));
bet=(Gamma1*eta1'+Gamma0*B*eta2')*100;

mu=rand(r,1);

X=floor(rand(n,p)*2)*100;
for i=1:n
    Y(i,:)=mu'+(bet*X(i,:)')'+mvnrnd(Mean,Sigma);
end
            
stat=ienv(X,Y,d);
subspace(stat.Gamma1,Gamma1)

alpha=0.01;
d=lrt_ienv(X,Y,alpha);
d=aic_ienv(X,Y);
d=bic_ienv(X,Y);

B=200;
bootse=bstrp_ienv(X,Y,B,d)


load irisf.mat

alpha=0.01;
u=bic_env(X,Y)
d=bic_ienv(X,Y)
stat=ienv(X,Y,d)
1-1./stat.ratio
%-----Test the scaled envelope model-----
r=10;
u=5;
p=5;
n=300;
sigma1=.5;
sigma2=sqrt(5);
Gamma=grams(rand(r,u));
Gamma0=grams(nulbasis(Gamma'));
bet=Gamma*randn(u,p)*2;
Sigmao=sigma1^2*Gamma*Gamma'+sigma2^2*Gamma0*Gamma0';
Mean=zeros(r,1);
X=randn(n,p)*5;

Z=zeros(n,r);
for i=1:n
    Z(i,:)=X(i,:)*bet'+mvnrnd(Mean,Sigmao);
end

seq=0:0.5:4.5;
D=diag(2.^seq);
Y=Z*D;
stat=senv(X,Y,u);
subspace(stat.Gamma,Gamma)
u=aic_senv(X,Y)
u=bic_senv(X,Y)

B=5;
bootse=bstrp_senv(X,Y,B,u)

%-----Test the heteroscedastic envelope model-----
sigma11=1;
sigma12=5;
sigma2=10;
p=2;
r=10;
u=2;
n=200;

Mean=zeros(r,1);
C1=grams(rand(r,r));
Gamma=C1(:,1:u);
Gamma0=C1(:,u+1:end);
Sigma1=sigma11^2*Gamma*Gamma'+sigma2^2*Gamma0*Gamma0';
Sigma2=sigma12^2*Gamma*Gamma'+sigma2^2*Gamma0*Gamma0';

eta=randn(u,1)*5;
bet=Gamma*eta;

mu=zeros(r,1);

Y=zeros(n,r);
X=ones(n,1);
X(1:n/2)=-1;

for i=1:n/2
    Y(i,:)=mu'+mvnrnd(Mean,Sigma1);
    Y(n/2+i,:)=mu'+bet'*2+mvnrnd(Mean,Sigma2);
end

B=10;
bootse=bstrp_henv(X,Y,B,u);

stat=henv(X,Y,0)
stat=henv(X,Y,r)
stat=henv(X,Y,u)

subspace(stat.Gamma,Gamma)
u=aic_henv(X,Y)
u=bic_henv(X,Y)
u=lrt_henv(X,Y,0.01)


% Waterstrider example
load waterstrider.mat
Y=log(Yu);
X=X([31:90 121:150],[2 3]);
X(61:90,:)=-1;
Y=Y([31:90 121:150],:);
k=3;
[n r]=size(Y);
u=lrt_henv(X,Y,0.01)
