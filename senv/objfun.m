function f=objfun(d,Gamma,dataParameter)

n=dataParameter.n;
sigY=dataParameter.sigY;
sigRes=dataParameter.sigRes;



    Lambda=diag([1 d]);
    eigtem = eig(Gamma'*Lambda*inv(sigY)*Lambda*Gamma);
    a= log(prod(eigtem(eigtem>0))); 
    eigtem2 = eig(Gamma'*inv(Lambda)*sigRes*inv(Lambda)*Gamma);
    b= log(prod(eigtem2(eigtem2>0)));
    
    f = n*a/2+n*b/2;