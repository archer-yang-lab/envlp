

function f = F4senv(R,dataParameter)

sigRes=dataParameter.sigRes;
sigY=dataParameter.sigY;
Lambda=dataParameter.Lambda;


eigtem=eig(R'*inv(Lambda)*sigRes*inv(Lambda)*R);
a=log(prod(eigtem(eigtem>0)));

eigtem0=eig(R'*Lambda*inv(sigY)*Lambda*R);
b=log(prod(eigtem0(eigtem0>0)));

f=a+b;