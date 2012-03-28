%% penv
% Fit the partial envelope model.

%% Usage
% stat=penv(X1,X2,Y,u)
%
% Input
%
% * X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.  The covariates can be univariate or 
% multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p1.
% * u: Dimension of the partial envelope. An integer between 0 and r.
%
% Output
% 
% stat: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * stat.beta1: The partial envelope estimator of $$\beta_1$, which is the
% regression coefficients for X1. An r by p1 matrix.
% * stat.beta2: The partial envelope estimator of $$\beta_2$, which is the
% regression coefficients for X2. An r by p2 matrix.
% * stat.Sigma: The partial envelope estimator of the error covariance 
% matrix.  An r by r matrix.
% * stat.Gamma: The orthogonal basis of the partial envelope subspace. An r by u
% semi-orthogonal matrix.
% * stat.Gamma0: The orthogonal basis of the complement of the partial envelope
% subspace.  An r by r-u semi-orthogonal matrix.
% * stat.eta: The coordinates of $$\beta_1$ with respect to Gamma. An u by
% p1 matrix.
% * stat.Omega: The coordinates of Sigma with respect to Gamma. An u by u
% matrix.
% * stat.Omega0: The coordinates of Sigma with respect to Gamma0. An r-u by r-u
% matrix.
% * stat.alpha: The estimated intercept in the partial envelope model.  An r by 1
% vector.
% * stat.l: The maximized log likelihood function.  A real number.
% * stat.asyPenv: Asymptotic standard error for elements in $$\beta$ under
% the partial envelope model.  An r by p1 matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1/sqrt(n).
% * stat.ratio: The asymptotic standard error ratio of the stanard multivariate 
% linear regression estimator over the partial envelope estimator, for each element 
% in $$\beta_1$.  An r by p1 matrix.
% * stat.np: The number of parameters in the envelope model.  A positive
% integer.

%% Description
% This function fits the partial envelope model to the responses Y and
% predictors X1 and X2, using the maximum likehood estimation.  When the dimension of the
% envelope is between 1 and r-1, we implemented the algorithm in Su and
% Cook (2011).  When the dimension is r, then the partial envelope model degenerates
% to the standard multivariate linear regression with Y as the responses and
% both X1 and X2 as predictors.  When the dimension is 0, X1 and Y are 
% uncorrelated, and the fitting is the standard multivariate linear
% regression with Y as the responses and X2 as the predictors.

%% References
% 
% * The codes is implemented based on the algorithm in Section 3.2 of Su
% and Cook (2012).
% * The Grassmann manifold optimization step calls the package sg_min 2.4.1
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).

%% Example
% 
% The following codes reconstruct the results in Su and Cook (2012).
% 
% load T7-7.dat
% Y=T7_7(:,1:4);
% X=T7_7(:,5:7);
% X1=X(:,3);
% X2=X(:,1:2);
% alpha=0.01;
% u=lrt_penv(X1,X2,Y,alpha)
% stat=penv(X1,X2,Y,u)
% stat.Omega
% eig(stat.Omega0)
% stat.ratio

function stat=penv(X1,X2,Y,u)

%---preparation---
[n p1]=size(X1);
p2=size(X2,2);
r=size(Y,2);

X1C=center(X1);
X2C=center(X2);
YC=center(Y);

SX1=cov(X1,1);
SX2=cov(X2,1);
SX12=X1C'*X2C/n;

QX2=eye(n)-X2C*inv(X2C'*X2C)*X2C';
R12=QX2*X1C;
RY2=QX2*YC;


% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the partial 
% envelope model reduces to the standard model, and it also should be fitted
% differently.


if u>0 && u<r


    %---Call env to compute most of the components in output---

    temp=env(R12,RY2,u);
    beta1=temp.beta;
    
    Gamma=temp.Gamma;
    Gamma0=temp.Gamma0;
    eta=temp.eta;
    Omega=temp.Omega;
    Omega0=temp.Omega0;
    Sigma1=Gamma*Omega*Gamma';
    Sigma=temp.Sigma;
    
    
    stat.beta1=beta1;
    stat.Gamma=Gamma;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;
    stat.Sigma=Sigma;
    stat.l=temp.l; 
    

    %---Compute the rest in output---
    beta2=(YC-X1C*beta1')'*X2C*inv(X2C'*X2C);
    alpha=mean(Y)'-beta1*mean(X1)'-beta2*mean(X2)';
    
    stat.beta2=beta2;    
    stat.alpha=alpha;
    stat.np=r+u*p1+r*p2+r*(r+1)/2;
    
    %---compute asymptotic variance and get the ratios---
    Sig1G2=SX1-SX12*inv(SX2)*SX12';
    asyFm=kron(inv(Sig1G2),Sigma);
    asyFm=reshape(sqrt(diag(asyFm)),r,p1);
    
    temp=kron(eta*Sig1G2*eta',inv(Omega0))+kron(Omega,inv(Omega0))+kron(inv(Omega),Omega0)-2*kron(eye(u),eye(r-u));
    asyPenv=kron(inv(Sig1G2),Sigma1)+kron(eta',Gamma0)*inv(temp)*kron(eta,Gamma0');
    asyPenv=reshape(sqrt(diag(asyPenv)),r,p1);
    
    stat.asyPenv=asyPenv;
    stat.ratio=asyFm./asyPenv;

    
elseif u==0
    
    [beta2 Sigma]=fit_OLS(X2,Y);
    eigtem=eig(Sigma);
    
    stat.beta1=zeros(r,p1);
    stat.beta2=beta2;
    stat.Sigma=Sigma;
    stat.eta=[];
    stat.Gamma=[];
    stat.Gamma0=eye(r);
    stat.Omega=[];
    stat.Omega0=Sigma;
    stat.alpha=mean(Y)'-beta2*mean(X2)';
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.asyPenv=[];
    stat.ratio=ones(r,p1);
    stat.np=r+u*p1+r*p2+r*(r+1)/2;
    

elseif u==r
    
    X=[X1 X2];
    [beta Sigma]=fit_OLS(X,Y);
    eigtem=eig(Sigma);
    sigX=cov(X,1);
    asyFm=kron(inv(sigX),Sigma);
    asyFm=reshape(sqrt(diag(asyFm(1:r*p1,1:r*p1))),r,p1);
    
    
    stat.beta1=beta(:,1:p1);
    stat.beta2=beta(:,p1+1:end);
    stat.Sigma=Sigma;
    stat.eta=beta(:,1:p1);
    stat.Gamma=eye(r);
    stat.Gamma0=[];
    stat.Omega=Sigma;
    stat.Omega0=[];
    stat.alpha=mean(Y)'-beta*mean(X)';
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.asyPenv=asyFm;
    stat.ratio=ones(r,p1);
    stat.np=r+u*p1+r*p2+r*(r+1)/2;
    
end
    