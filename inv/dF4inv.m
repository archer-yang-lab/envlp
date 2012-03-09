function df = dF4inv(R,dataParameter)

[r d]=size(R);

sigRes=dataParameter.sigRes;
sigY=dataParameter.sigY;
sigFit=sigY-sigRes;
p=dataParameter.p;

if d==p
    
    a=sigRes*R*inv(R'*sigRes*R);
    temp=inv(sigY);
    b=temp*R*inv(R'*temp*R);
    df=a+b;
    
else
    
    temp=inv(sigRes);
    a=2*sigRes*R*inv(R'*sigRes*R)+2*temp*R*inv(R'*temp*R);

    temp1=inv(R0'*sigRes*R0);
    temp2=R0'*sigFit*R0;
    R0=grams(nulbasis(R'));
    dzdg0=kron(eye(r-d),temp1*R0'*sigFit)+Kpd(kron(temp1,R0'*sigFit),r-d,r-d)-kron(temp2*temp1,temp1*R0'*sigRes)-Kpd(kron(temp1,temp2*temp1*R0'*sigRes),r-d,r-d);
    dg0dg1=-Kpd_right(kron(R0',R),r,d);
    
    [V0 D0]=eig(temp1);
    MultiPlier=V*diag(sqrt(diag(D)))*V';
    [Vtmp Dtmp]=eig(MultiPlier*temp2*MultiPlier);
    [Ds ind]=sort(diag(Dtmp),'descend');
    b=zeros(1,d*r);

    
    for i=(p-d+1):(r-d)
        V=inv(MultiPlier)*Vtmp(:,i);
        U=MultiPlier*Vtmp(:,i);
        b=b+1/(1+Ds(i))*reshape(V(:,i)*U(:,i)',1,r^2);
    end
    
    b=b*dzdg0*dg0dg1;
    
    df=a+reshape(b,r,d);
end