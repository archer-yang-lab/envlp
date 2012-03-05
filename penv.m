function stat=penv(X1,X2,Y,u)

%---preparation---
[n p1]=size(X1);
p2=size(X2,2);
r=size(Y,2);

X1C=center(X1);
X2C=center(X2);
YC=center(Y);

SX1=cov(X1,1);
SX2=cov(X2,1);
SX12=X1C'*X2C/n;

QX2=eye(n)-X2C*inv(X2C'*X2C)*X2C';
R12=QX2*X1C;
RY2=QX2*YC;


% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the partial 
% envelope model reduces to the standard model, and it also should be fitted
% differently.


if u>0 && u<r


    %---Call env to compute most of the components in output---

    temp=env(R12,RY2,u);
    beta1=temp.beta;
    
    Gamma=temp.Gamma;
    Gamma0=temp.Gamma0;
    eta=temp.eta;
    Omega=temp.Omega;
    Omega0=temp.Omega0;
    Sigma1=Gamma*Omega*Gamma';
    Sigma=temp.Sigma;
    
    
    stat.beta1=beta1;
    stat.Gamma=Gamma;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;
    stat.Sigma=Sigma;
    stat.l=temp.l; 
    

    %---Compute the rest in output---
    beta2=(YC-X1C*beta1')'*X2C*inv(X2C'*X2C);
    alpha=mean(Y)'-beta1*mean(X1)'-beta2*mean(X2)';
    
    stat.beta2=beta2;    
    stat.alpha=alpha;
    stat.np=r+u*p1+r*p2+r*(r+1)/2;
    
    %---compute asymptotic variance and get the ratios---
    Sig1G2=SX1-SX12*inv(SX2)*SX12';
    asyfm=kron(inv(Sig1G2),Sigma);
    temp=kron(eta*Sig1G2*eta',inv(Omega0))+kron(Omega,inv(Omega0))+kron(inv(Omega),Omega0)-2*kron(eye(u),eye(r-u));
    asypm=kron(inv(Sig1G2),Sigma1)+kron(eta',Gamma0)*inv(temp)*kron(eta,Gamma0');
    stat.ratio=reshape(sqrt(diag(asyfm)./diag(asypm)),r,p1);

    
elseif u==0
    
    [beta2 Sigma]=fit_OLS(X2,Y);
    eigtem=eig(Sigma);
    
    stat.beta1=zeros(r,p1);
    stat.beta2=beta2;
    stat.Sigma=Sigma;
    stat.eta=[];
    stat.Gamma=[];
    stat.Gamma0=eye(r);
    stat.Omega=[];
    stat.Omega0=Sigma;
    stat.alpha=mean(Y)'-beta2*mean(X2)';
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.ratio=ones(r,p1);
    stat.np=r+u*p1+r*p2+r*(r+1)/2;
    

elseif u==r
    
    X=[X1 X2];
    [beta Sigma]=fit_OLS(X,Y);
    eigtem=eig(Sigma);
    
    stat.beta1=beta(:,1:p1);
    stat.beta2=beta(:,p1+1:end);
    stat.Sigma=Sigma;
    stat.eta=beta(:,1:p1);
    stat.Gamma=eye(r);
    stat.Gamma0=[];
    stat.Omega=Sigma;
    stat.Omega0=[];
    stat.alpha=mean(Y)'-beta*mean(X)';
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.ratio=ones(r,p1);
    stat.np=r+u*p1+r*p2+r*(r+1)/2;
    
end
    