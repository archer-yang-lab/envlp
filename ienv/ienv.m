function stat=ienv(X,Y,u)


dataParameter=make_parameter(X,Y);
n=dataParameter.n;
p=dataParameter.p;
r=dataParameter.r;
XC=dataParameter.XC;
YC=dataParameter.YC;
mX=dataParameter.mX;
mY=dataParameter.mY;
sigX=dataParameter.sigX;
sigY=dataParameter.sigY;
sigRes=dataParameter.sigRes;
betaOLS=dataParameter.betaOLS;

if r<=p
    
    error('When the number of responses is less than the number of predictors, the inner envelope model cannot be applied.');
    
end

if u>p
    
    error('The dimension of the inner envelope subspace cannot be greater than p.');

elseif u==p
    
    temp=env(X,Y,p);
    stat.beta=temp.beta;
    stat.Sigma=temp.Sigma;
    stat.Gamma1=temp.Gamma;
    stat.Gamma0=temp.Gamma0;
    stat.B=[];
    stat.eta1=temp.eta';
    stat.eta2=[];
    stat.Omega1=temp.Omega;
    stat.Omega0=temp.Omega0;
    stat.alpha=temp.alpha;
    stat.l=temp.l;
    stat.np=temp.np;
    stat.ratio=temp.ratio;
    
elseif u==0
    
    eigtem=eig(sigY);
    stat.beta=zeros(r,p);
    stat.Sigma=sigY;
    stat.Gamma1=[];
    stat.Gamma0=eye(r);
    stat.B=eye(r);
    stat.eta1=[];
    stat.eta2=[];
    stat.Omega=[];
    stat.Omega0=sigY;
    stat.alpha=mY;
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.ratio=ones(r,p);
    stat.np=r+u*p+r*(r+1)/2;    

else 
    
    sigFit=sigY-sigRes;
    eigtem=eig(sigRes);
    
    F = make_F(@F4ienv,dataParameter);
    dF = make_dF(@dFienv,dataParameter);
    
    init=get_Init(F,X,Y,u,dataParameter);
    [l Gamma1]=sg_min(F,dF,init,'prcg','quiet');
    
    
    
    Gamma0=grams(nulbasis(Gamma1'));
    eta1=(Gamma1'*betaOLS)';
    Omega1=(YC*Gamma1-XC*eta1)'*(YC*Gamma1-XC*eta1)/n;
    
    [Vtmp Dtmp]=eig(inv(Gamma0'*sigRes*Gamma0));
    temp1=Vtmp*diag(sqrt(diag(Dtmp)))*Vtmp';
    temp2=Gamma0'*sigFit*Gamma0;
    temp3=Vtmp*diag(diag(Dtmp).^(-0.5))*Vtmp';
    [Vt Lambdat]=eig(temp1*temp2*temp1);
    [Lambda ind]=sort(diag(Lambdat),'descend');
    V=Vt(:,ind);
    K=diag([zeros(1,p-u),Lambda(p-u+1:end)]);
    Omega0=Gamma0'*sigRes*Gamma0+temp3*V*K*V'*temp3;
    
    [Vtmp Dtmp]=eig(inv(Omega0));
    Omega0inv12=Vtmp*diag(sqrt(diag(Dtmp)))*Vtmp;
    [Vtmp Dtmp]=eig(Omega0inv12*Gamma0'*sigFit*Gamma0*Omega0inv12);
    [Ds ind]=sort(Dtmp,'descend');
    B=grams(inv(Omega0inv12)*Vtmp(:,ind(1:p-u)));
    eta2=(inv(B'*inv(Omega0)*B)*B'*inv(Omega0)*Gamma0'*betaOLS)';
    beta=Gamma1*eta1'+Gamma0*B*eta2';
    Sigma=Gamma1*Omega1*Gamma1'+Gamma0*Omega0*Gamma0';
    alpha=mY-beta*mX;
    
    stat.beta=beta;
    stat.Sigma=Sigma;
    stat.Gamma1=Gamma1;
    stat.Gamma0=Gamma0;
    stat.B=B;
    stat.eta1=eta1;
    stat.eta2=eta2;
    stat.Omega1=Omega1;
    stat.Omega0=Omega0;
    stat.alpha=alpha;
    stat.np=p^2+(p-u)*(r-p)+r*(r+1)/2;
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)))-n/2*l;
    
    
    %-----Compute asymptotic variance for inner envelope model-----
    %-----Standard Model-----
    asyfm=kron(inv(sigX),Sigma);
    
    %-----Inner Envelope Model-----
    B0=grams(nulbasis(B'));

    H=zeros(p*r+(r+1)*r/2,p*p+(p-u)*(r-p)+u*r+u*(u+1)/2+(r-u)*(r-u+1)/2);

    H(1:p*r,1:p*u)=kron(eye(p),Gamma1);
    H(1:p*r,p*u+1:p*p)=kron(eye(p),Gamma0*B);
    H(1:p*r,p*p+1:p*p+(p-u)*(r-p))=kron(eta2,Gamma0*B0);
    H(1:p*r,p*p+(p-u)*(r-p)+1:p*p+(p-u)*(r-p)+r*u)=kron(eta1,eye(r))-Kpd_right(kron(eta2*B'*Gamma0',Gamma1),r,u);
    H(p*r+1:end,p*p+(p-u)*(r-p)+1:p*p+(p-u)*(r-p)+r*u)=2*Contr(r)*(kron(Gamma1*Omega1,eye(r))-kron(Gamma1,Gamma0*Omega0*Gamma0'));
    H(p*r+1:end,p*p+(p-u)*(r-p)+r*u+1:p*p+(p-u)*(r-p)+r*u+u*(u+1)/2)=Contr(r)*kron(Gamma1,Gamma1)*Expan(u);
    H(p*r+1:end,p*p+(p-u)*(r-p)+r*u+u*(u+1)/2+1:end)=Contr(r)*kron(Gamma0,Gamma0)*Expan(r-u);

    insigma=inv(Sigma);
    J=zeros(p*r+(r+1)*r/2,p*r+(r+1)*r/2);
    J(1:p*r,1:p*r)=kron(Sx,insigma);
    J(p*r+1:end,p*r+1:end)=Expan(r)'*kron(insigma,insigma)*Expan(r)/2;
    asyv=H*pinv(H'*J*H)*H';
    asyim=asyv(r*p,r*p);
    stat.ratio=reshape(sqrt(diag(asyfm)./diag(asyim)),r,p);
    
end


