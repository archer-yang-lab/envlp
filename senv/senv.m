%% senv
% Fit the scaled envelope model.

%% Usage
% stat=senv(X,Y,u)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p.
% * u: Dimension of the envelope. An integer between 0 and r.
%
% Output
% 
% stat: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * stat.beta: The scaled envelope estimator of the regression coefficients
% $$\beta$. An r by p matrix.
% * stat.Sigma: The scaled envelope estimator of the error covariance
% matrix.  An r by r matrix.
% * stat.Lambda: The matrix of estimated scales. An r by r diagonal matrix
% with the first diagonal element equal to 1 and other diagonal elements
% being positive.
% * stat.Gamma: The orthogonal basis of the envelope subspace. An r by u
% semi-orthogonal matrix.
% * stat.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An r by r-u semi-orthogonal matrix.
% * stat.eta: The coordinates of $$\beta$ with respect to Gamma. An u by p
% matrix.
% * stat.Omega: The coordinates of Sigma with respect to Gamma. An u by u
% matrix.
% * stat.Omega0: The coordinates of Sigma with respect to Gamma0. An r-u by r-u
% matrix.
% * stat.alpha: The estimated intercept in the scaled envelope model.  An r
% by 1 vector.
% * stat.l: The maximized log likelihood function.  A real number.
% * stat.asySenv: Asymptotic standard error for elements in $$\beta$ under
% the scaled envelope model.  An r by p matrix.
% * stat.ratio: The asymptotic standard error ratio of the standard
% multivariate linear regression estimator over the scaled envelope
% estimator, for each element in $$\beta$.  An r by p matrix.
% * stat.np: The number of parameters in the scaled envelope model.  A
% positive integer.

%% Description
% This function fits the scaled envelope model to the responses and predictors,
% using the maximum likehood estimation.  When the dimension of the
% envelope is between 1 and r-1, we implemented the algorithm in Cook and
% Su (2012).  When the dimension is r, then the scaled envelope model 
% degenerates to the standard multivariate linear regression.  When the
% dimension is 0, it means that X and Y are uncorrelated, and the fitting
% is different.

%% References
% 
% * The codes is implemented based on the algorithm in Section 4.1 of Cook 
% and Su (2012).
% * The Grassmann manifold optimization step calls the package sg_min 2.4.1
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).


function stat=senv(X,Y,u)


dataParameter=make_parameter(X,Y);
n=dataParameter.n;
p=dataParameter.p;
r=dataParameter.r;
mX=dataParameter.mX;
mY=dataParameter.mY;
sigX=dataParameter.sigX;
sigY=dataParameter.sigY;
sigRes=dataParameter.sigRes;
betaOLS=dataParameter.betaOLS;

if u==0
    
    eigtem=eig(sigY);
    
    stat.beta=zeros(r,p);
    stat.Sigma=sigY;
    stat.Lambda=eye(r);
    stat.Gamma=[];
    stat.Gamma0=eye(r);
    stat.eta=[];
    stat.Omega=[];
    stat.Omega0=sigY;
    stat.alpha=mY;
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.asySenv=[];
    stat.ratio=ones(r,p);
    stat.np=r+u*p+r*(r+1)/2;  
    
elseif u==r
    
    asyFm=kron(inv(sigX),sigRes);
    asyFm=reshape(sqrt(diag(asyFm)),r,p);
    
    eigtem=eig(sigRes);
    
    stat.beta=betaOLS;
    stat.Sigma=sigRes;
    stat.Lambda=eye(r);
    stat.Gamma=eye(r);
    stat.Gamma0=[];
    stat.eta=betaOLS;
    stat.Omega=sigRes;
    stat.Omega0=[];
    stat.alpha=mY-betaOLS*mX;
    stat.l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.asySenv=asyFm;
    stat.ratio=ones(r,p);
    stat.np=r+u*p+r*(r+1)/2;
    
else

    ite=1000;
    epsilon=1e-9; 
    l2=zeros(1,ite);
    
    init=env(X,Y,u);
    Gamma=init.Gamma;
    d=ones(1,r-1);

    
    
    for i=1:ite

        Lambda=diag([1 d]);
        dataParameter.Lambda=Lambda;  
        
        [d l2(i)]=fminsearch(@(d) objfun(d,Gamma,dataParameter), d);
        
        if (i>1) && (abs(l2(i)-l2(i-1))<epsilon*abs(l2(i)))
            break;
        end
        
        F = make_F(@F4senv,dataParameter);
        dF = make_dF(@dF4senv,dataParameter); 
        [l Gamma]=sg_min(F,dF,Gamma,'prcg','quiet');
    
    end
    
    Lambda=diag([1 d]);
    Gamma0=grams(nulbasis(Gamma'));
    eta=Gamma'*inv(Lambda)*betaOLS;
    beta=Lambda*Gamma*eta;
    Omega=Gamma'*inv(Lambda)*sigRes*inv(Lambda)*Gamma;
    Omega0=Gamma0'*inv(Lambda)*sigY*inv(Lambda)*Gamma0;    
    Sigma=Lambda*(Gamma*Omega*Gamma'+Gamma0*Omega0*Gamma0')*Lambda;

    eigtem=eig(sigY);
    
    stat.beta=beta;
    stat.Sigma=Sigma;
    stat.Lambda=Lambda;
    stat.Gamma=Gamma;
    stat.Gamma0=Gamma0;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;    
    stat.alpha=mY-beta*mX;
    stat.np=init.np+r-1;
    stat.l=-n*r/2*(1+log(2*pi))-n/2*(l+log(prod(eigtem(eigtem>0))));
    
   %---compute asymptotic variance and get the ratios---
    asyFm=kron(inv(sigX),Sigma);
    asyFm=reshape(sqrt(diag(asyFm)),r,p);
    insigma=inv(Sigma);
        
    J=zeros(p*r+(r+1)*r/2,p*r+(r+1)*r/2);
    J(1:p*r,1:p*r)=kron(sigX,insigma);
    J(p*r+1:end,p*r+1:end)=Expan(r)'*kron(insigma,insigma)*Expan(r)/2;
    
    H=zeros(p*r+(r+1)*r/2,r-1+p*u+r*(r+1)/2);

    H(1:p*r,1:r-1)=kron(eta'*Gamma',eye(r))*Lmatrix(r);
    H(1:p*r,r:r-1+p*u)=kron(eye(p),Lambda*Gamma);
    H(1:p*r,r+p*u:r-1+u*(r-u+p))=kron(eta',Lambda*Gamma0);

    H(p*r+1:end,1:r-1)=Contr(r)*(kron(Lambda*Gamma*Omega*Gamma',eye(r))+kron(eye(r),Lambda*Gamma*Omega*Gamma'))*Lmatrix(r)+Contr(r)*(kron(Lambda*Gamma0*Omega0*Gamma0',eye(r))+kron(eye(r),Lambda*Gamma0*Omega0*Gamma0'))*Lmatrix(r);
    H(p*r+1:end,r+p*u:r-1+u*(r-u+p))=2*Contr(r)*(kron(Lambda*Gamma*Omega,Lambda*Gamma0)-kron(Lambda*Gamma,Lambda*Gamma0*Omega0));
    H(p*r+1:end,r+u*(r-u+p):r-1+u*(r-u+p)+u*(u+1)/2)=Contr(r)*kron(Lambda*Gamma,Lambda*Gamma)*Expan(u);
    H(p*r+1:end,r+u*(r-u+p)+u*(u+1)/2:end)=Contr(r)*kron(Lambda*Gamma0,Lambda*Gamma0)*Expan(r-u);

    asyvar=H*inv(H'*J*H)*H'; 
    asySenv=reshape(sqrt(diag(asyvar(1:r*p,1:r*p))),r,p);
    stat.asySenv=asySenv;
    stat.ratio=asyFm./asySenv;
    
end