function [beta_OLS Sigma_OLS]=fit_OLS(X,Y)

n=length(X);
XC=center(X);
YC=center(Y);
PX=XC*inv(XC'*XC)*XC';
beta_OLS=YC'*XC*inv(XC'*XC);
Sigma_OLS=YC'*(eye(n)-PX)*YC/n;