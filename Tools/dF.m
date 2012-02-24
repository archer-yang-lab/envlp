function df = dF(R)

global sigY;
global sigres;

a=sigres*R*inv(R'*sigres*R);

temp=inv(sigY);

b=temp*R*inv(R'*temp*R);

df=a+b;