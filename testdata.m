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
bet=G*ita;
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

u=lrt_env(Y,X,alpha)
u=bic_env(Y,X)
u=aic_env(Y,X)



load wheatprotein.txt
X=wheatprotein(:,8);
Y=wheatprotein(:,1:6);
u=lrt_env(Y,X,alpha)
u=bic_env(Y,X)
u=aic_env(Y,X)



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
u=aic_penv(X1,X2,Y)
stat=penv(X1,X2,Y,1)

X1=X(:,2);
X2=X(:,[1 3]);
