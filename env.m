function [beta Sigma Gamma Gamma0 eta Omega Omega0 alpha l ratio]=env(X,Y,u)


% To Yi: 1) Maybe we can return a list that includes all the parameters, 
% which makes the interface cleaner. But I do not know how to do it. 2) We
% need to do check something, e.g., if Y is discrete, the model cannot
% handle that, but I do not know how to check that either.  There are some
% other checks: u must be an interger between 0 and r, X and Y must have
% the same length.


global sigY;
global sigres;

%---preparation---
[n p]=size(X);
r=size(Y,2);
XC=center(X);
YC=center(Y);

sigX=cov(X,1);
sigY=cov(Y,1);

[beta_OLS sigres]=fit_OLS(X,Y);
eigtem=eig(sigY);


% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the envelope
% model reduces to the standard model, and it also should be fitted
% differently.


if u>0 && u<r


    %---Compute \Gamma using sg_min---

    init=startv(X,Y,u);
    [l Gamma]=sg_min(init,'prcg','quiet');


    %---Compute the rest of the parameters based on \Gamma---
    Gamma0=grams(nulbasis(Gamma'));
    alpha=mean(Y)';
    beta=Gamma*Gamma'*beta_OLS;
    eta=Gamma'*beta;
    Omega=Gamma'*sigres*Gamma;
    Omega0=Gamma0'*sigY*Gamma0;
    Sigma1=Gamma*Omega*Gamma';
    Sigma2=Gamma0*Omega0*Gamma0';
    Sigma=Sigma1+Sigma2;
    l=-n*r/2*(1+log(2*pi))-n/2*(l+log(prod(eigtem(eigtem>0))));

    %---compute asymptotic variance and get the ratios---
    asyfm=kron(inv(cov(X,1)),Sigma);
    temp=kron(eta*sigX*eta',inv(Omega0))+kron(Omega,inv(Omega0))+kron(inv(Omega),Omega0)-2*kron(eye(u),eye(r-u));
    asyem=kron(inv(sigX),Sigma1)+kron(eta',Gamma0)*inv(temp)*kron(eta,Gamma0');
    ratio=reshape(sqrt(diag(asyfm)./diag(asyem)),r,p);

    
    
elseif u==0
    
    
    
    beta=zeros(r,p);
    Gamma=[];
    eta=[];
    Omega=[];
    Gamma0=eye(r);
    Sigma=sigY;
    Omega0=sigY;
    alpha=mean(Y)';
    l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    ratio=ones(r,p);
    
    

elseif u==r
    
    
    beta=beta_OLS;
    eta=beta;
    Sigma=sigres;
    Gamma=eye(r);
    Gamma0=[];
    Omega=sigres;
    Omega0=[];
    alpha=mean(Y)';
    eigtem=eig(sigres);
    l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    ratio=ones(r,p);
    
end
    
    
    