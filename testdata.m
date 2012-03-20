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



load wheatprotein.txt
X=wheatprotein(:,8);
Y=wheatprotein(:,1:6);
alpha=0.01;
u=lrt_env(Y,X,alpha)
u=bic_env(Y,X)
u=aic_env(Y,X)
stat=env(X,Y,u);

%-----Test partial envelope model-----

load T7-7.dat
Y=T7_7(:,1:4);
X=T7_7(:,5:7);
X1=X(:,3);
X2=X(:,1:2);
u=lrt_env(Y,X,0.01);
stat=env(X,Y,u);
stat=penv(X1,X2,Y,1)

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