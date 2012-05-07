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
ModelOutput=env(X,Y,u);
eig(ModelOutput.Omega)
eig(ModelOutput.Omega0)
norm(ModelOutput.beta-bet)
ModelOutput=fit_OLS(X,Y);
beta_OLS=ModelOutput.betaOLS;
sigres=ModelOutput.SigmaOLS;
norm(beta_OLS-bet)
subspace(ModelOutput.Gamma,G)

alpha=0.01;
u=lrt_env(X,Y,alpha)
u=bic_env(X,Y)
u=aic_env(X,Y)



bootse=bstrp_OLS(X,Y,B)
bootse=bstrp_env(X,Y,B,u)

load wheatprotein.txt
X=wheatprotein(:,8);
Y=wheatprotein(:,1:6);
bootse=bstrp_OLS(X,Y,200)


alpha=0.01;
u=lrt_env(X,Y,alpha)
ModelOutput=env(X,Y,u);
ModelOutput.Omega
eig(ModelOutput.Omega0)
ModelOutput.ratio

u=bic_env(X,Y)
u=aic_env(X,Y)


%-----Test partial envelope model-----

load T7-7.dat
Y=T7_7(:,1:4);
X=T7_7(:,5:7);
X1=X(:,3);
X2=X(:,1:2);
alpha=0.01;
u=lrt_penv(X1,X2,Y,alpha)
u=lrt_env(X,Y,0.01);
ModelOutput=env(X,Y,u);
ModelOutput=penv(X1,X2,Y,1)


B=200;
bootse=bstrp_penv(X1,X2,Y,B,u)

X1=X(:,3);
X2=X(:,1:2);
u=lrt_penv(X1,X2,Y,0.01)
ModelOutput=penv(X1,X2,Y,u)

X1=X(:,2);
X2=X(:,[1 3]);
u=lrt_penv(X1,X2,Y,0.01)
ModelOutput=penv(X1,X2,Y,u)

X1=X(:,1);
X2=X(:,2:3);
u=lrt_penv(X1,X2,Y,0.01)
ModelOutput=penv(X1,X2,Y,u)


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
            
ModelOutput=ienv(X,Y,d);
subspace(ModelOutput.Gamma1,Gamma1)

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
ModelOutput=ienv(X,Y,d)
1-1./ModelOutput.ratio
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
ModelOutput=senv(X,Y,u);
subspace(ModelOutput.Gamma,Gamma)
u=aic_senv(X,Y)
u=bic_senv(X,Y)

B=5;
bootse=bstrp_senv(X,Y,B,u)


load('T9-12.txt')
Y=T9_12(:,4:7);
X=T9_12(:,1:3);
u=bic_env(X,Y)
ModelOutput=env(X,Y,u);
1-1./ModelOutput.ratio
u=bic_senv(X,Y)
ModelOutput=senv(X,Y,u);
ModelOutput.Lambda
1-1./ModelOutput.ratio
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

ModelOutput=henv(X,Y,0)
ModelOutput=henv(X,Y,r)
ModelOutput=henv(X,Y,u)

subspace(ModelOutput.Gamma,Gamma)
u=aic_henv(X,Y)
u=bic_henv(X,Y)
u=lrt_henv(X,Y,0.01)


% Waterstrider example
load waterstrider.mat
[n r]=size(Y);
u=lrt_henv(X,Y,0.01)
ModelOutput=henv(X,Y,u)
ModelOutput.ratio


% Test xenv
load wheatprotein.txt
X=wheatprotein(:,1:6);
Y=wheatprotein(:,7);
ModelOutput=xenv(X,Y,0);

p=size(X,2);
ModelOutput=xenv(X,Y,p);

% When u=p, the envelope model reduces to the ordinary least squares
% regression

temp=fit_OLS(X,Y);
temp.SigmaOLS
ModelOutput.sigYcX
temp.betaOLS'
ModelOutput.beta

ModelOutput=xenv(X,Y,1);
u=aic_xenv(X,Y)
u=bic_xenv(X,Y)