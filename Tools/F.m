function f = F(R)

global sigY;
global sigres;

eigtem=eig(R'*sigres*R);
a=log(prod(eigtem(eigtem>0)));

eigtem0=eig(R'*inv(sigY)*R);
b=log(prod(eigtem(eigtem>0)));

f=a+b;