function dataParameter=make_parameter(X,Y)

[n p]=size(X);
r=size(Y,2);

XC=center(X);
YC=center(Y);
[betaOLS sigRes]=fit_OLS(X,Y);

dataParameter.n=n;
dataParameter.p=p;
dataParameter.r=r;
dataParameter.XC=XC;
dataParameter.YC=YC;
dataParameter.mX=mean(X)';
dataParameter.mY=mean(Y)';
dataParameter.sigX=cov(X,1);
dataParameter.sigY=cov(Y,1);
dataParameter.sigRes=sigRes;
dataParameter.betaOLS=betaOLS;
