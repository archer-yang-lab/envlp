n=100;
r=10;
u=1;
p=2;

X=rand(n,p);
GG=rand(r,r);
GG=grams(GG);
G=GG(:,1:u);
G0=GG(:,u+1:end);
sigma=1;
sigma0=5;
eta=rand(u,p);
bet=G*eta;
Sigma=G*G'*sigma^2+G0*G0'*sigma0^2;
epsil=mvnrnd(zeros(1,r),Sigma,n);
Y=X*bet'+epsil;

[beta Sigma Gamma Gamma0 eta Omega Omega0 alpha l ratio]=env(X,Y,u)