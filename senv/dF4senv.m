
function df = dF4senv(R,dataParameter)

sigRes=dataParameter.sigRes;
sigY=dataParameter.sigY;
Lambda=dataParameter.Lambda;

a=2*inv(Lambda)*sigRes*inv(Lambda)*R*inv(R'*inv(Lambda)*sigRes*inv(Lambda)*R);

temp=inv(sigY);

b=2*Lambda*temp*Lambda*R*inv(R'*Lambda*temp*Lambda*R);

df=a+b;