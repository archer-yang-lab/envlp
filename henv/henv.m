%% henv
% Fit the heteroscedastic envelope model.

%% Usage
% stat=henv(X,Y,u)
%
% Input
%
% * X: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
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
% * stat.mu: The heteroscedastic envelope estimator of the grand mean. A r
% by 1 vector.
% * stat.mug: The heteroscedastic envelope estimator of the group mean. A r
% by p matrix, the ith column of the matrix contains the mean for the ith
% group.
% * stat.Yfit: A n by r matrix, the ith row gives the group mean of the
% group that the ith observation belongs to.  As X is just a group
% indicator, and is not ordinal, stat.mug alone does not tell which
% group corresponds to which group mean.
% * stat.Gamma: The orthogonal basis of the envelope subspace. An r by u
% semi-orthogonal matrix.
% * stat.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An r by r-u semi-orthogonal matrix.
% * stat.beta: The heteroscedastic envelope estimator of the group mean
% effect. An r by p matrix, the ith column of the matrix contains the
% main effect for the ith group.
% * stat.Sigma: The heteroscedastic envelope estimator of the error
% covariance matrix.  A three dimensional matrix with dimension r, r and p,
% stat.Sigma(:,:,i) contains the estimated covariance matrix for the ith
% group.
% * stat.eta: The coordinates of $$\beta$ with respect to Gamma. An u by p
% matrix, the ith column contains the coordinates of the main effect of the
% ith group with respect to Gamma.
% * stat.Omega: The coordinates of Sigma with respect to Gamma. An u by u
% by p matrix, stat.Omega(:,:,i) contains the coordinates of the covariance
% matrix of the ith group with respect to Gamma.
% * stat.Omega0: The coordinates of Sigma with respect to Gamma0. An r-u by r-u
% matrix.
% * stat.l: The maximized log likelihood function.  A real number.
% * stat.np: The number of parameters in the heteroscedastic envelope
% model.  A positive integer.
% * stat.asyHenv: The asymptotic standard errors for elements in $$beta$
% under the heteroscedastic envelope model. An r by p matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1/sqrt(n).
% * stat.ratio: The asymptotic standard error ratio of the standard multivariate 
% linear regression estimator over the heteroscedastic envelope estimator.
% An r by p matrix, the (i, j)th element in stat.ratio is the elementwise standard
% error ratio for the ith element in the jth group mean effect.


%% Description
% This function fits the heteroscedatic envelope model to the responses and predictors,
% using the maximum likehood estimation.  When the dimension of the
% envelope is between 1 and r-1, we implemented the algorithm in Su and Cook (2012).
% When the dimension is r, then the envelope model degenerates
% to the standard multivariate linear model for comparing group means.  When the dimension is 0,
% it means there is not any group effect, and the fitting is different.

%% References
% 
% * The codes is implemented based on the algorithm in Section 2.2 of Su
% and Cook (2012).
% * The Grassmann manifold optimization step calls the package sg_min 2.4.1
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).


function stat=henv(X,Y,u)

dataParameter=make_parameter4henv(X,Y);

p=dataParameter.p;
r=dataParameter.r;
n=dataParameter.n;
ng=dataParameter.ng;
ncum=dataParameter.ncum;
mY=dataParameter.mY;
mYg=dataParameter.mYg;
sigRes=dataParameter.sigRes;
sigY=dataParameter.sigY;
ind=dataParameter.ind;

if u==0
    
    Gamma=[];
    Gamma0=eye(r);
    Sigma=sigY;
    mu=mY;
    mug=mY*ones(1,p);
    beta=zeros(r,p);
    eta=[];
    Omega=[];
    Omega0=sigY;
    eigtem = eig(sigY);
    l=-n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    
    stat.mu=mu;
    stat.mug=mug;
    stat.Yfit=zeros(n,r);
    stat.Gamma=Gamma;
    stat.Gamma0=Gamma0;
    stat.beta=beta;
    stat.Sigma=Sigma;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;
    stat.l=l;
    stat.np=(r-u)+u*(r-u+p)+p*u*(u+1)/2+(r-u)*(r-u+1)/2;
    stat.asyHenv=[];
    stat.ratio=ones(r,p);    
    
    
elseif u==r
    
    Gamma=eye(r);
    Gamma0=[];
    Sigma=sigRes;
    mu=mY;
    mug=mYg;
    beta=mug-mu*ones(1,p);
    eta=beta;
    Omega=sigRes;
    Omega0=[];
    l=-n*r/2*(1+log(2*pi));
    for i=1:p
        eigtem = eig(sigRes(:,:,i));
        l=l-ng(i)/2*log(prod(eigtem(eigtem>0)));
    end
    Yfit=zeros(n,r);
    for i=1:p
        if i>1
            Yfit(ind(ncum(i-1)+1:ncum(i)),:)=ones(ng(i),1)*mug(:,i)';
        else
            Yfit(ind(1:ncum(1)),:)=ones(ng(1),1)*mug(:,1)';
        end
    end
    
    fracN=ng/n;
    J=zeros(p*r+p*r*(r+1)/2);
    for i=1:p-1
        for j=1:p-1

            J((i-1)*r+1:i*r,(j-1)*r+1:j*r)=fracN(j)*fracN(i)/fracN(p)*inv(Sigma(:,:,p));
        end
        J((i-1)*r+1:i*r,(i-1)*r+1:i*r)=fracN(i)*inv(Sigma(:,:,i))+fracN(i)^2/fracN(p)*inv(Sigma(:,:,p));
    end
    for i=1:p
        J(r*(p-1)+1+(i-1)*r*(r+1)/2:r*(p-1)+i*r*(r+1)/2,r*(p-1)+1+(i-1)*r*(r+1)/2:r*(p-1)+i*r*(r+1)/2)=0.5*fracN(i)*Expan(r)'*kron(inv(Sigma(:,:,i)),inv(Sigma(:,:,i)))*Expan(r);
    end
    J(r+1:end,r+1:end)=J(1:(p-1)*r+p*r*(r+1)/2,1:(p-1)*r+p*r*(r+1)/2);
    J(1:r,:)=0;
    J(r+1:end,1:r)=0;
    for i=1:p
        J(1:r,1:r)=J(1:r,1:r)+fracN(i)*inv(Sigma(:,:,i));
    end
    for i=1:p-1
        J(1:r,i*r+1:(i+1)*r)=fracN(i)*(inv(Sigma(:,:,p))-inv(Sigma(:,:,i)));
        J(i*r+1:(i+1)*r,1:r)=J(1:r,i*r+1:(i+1)*r);
    end
    temp=inv(J);
    asyFm=sqrt(diag(temp(1:r*p,1:r*p)));
    
    stat.mu=mu;
    stat.mug=mug;
    stat.Yfit=Yfit;
    stat.Gamma=Gamma;
    stat.Gamma0=Gamma0;
    stat.beta=beta;
    stat.Sigma=Sigma;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;
    stat.l=l;
    stat.np=(r-u)+u*(r-u+p)+p*u*(u+1)/2+(r-u)*(r-u+1)/2;
    stat.asyHenv=reshape(asyFm,r,p);
    stat.ratio=ones(r,p);
    
else
    
    mu=mY;
    eigtem = eig(sigY);
 
    F = make_F(@F4henv,dataParameter);
    dF = make_dF(@dF4henv,dataParameter);
    
    init=get_Init4henv(F,X,Y,u,dataParameter);
    [l Gamma]=sg_min(F,dF,init,'prcg','quiet');

    Gamma0=grams(nulbasis(Gamma'));
    Omega0=Gamma0'*sigY*Gamma0;
    eta=Gamma'*(mYg-mu*ones(1,p));
    beta=Gamma*eta;
    mug=mu*ones(1,p)+beta;
    Omega=zeros(u,u,p);
    Sigma=zeros(r,r,p);
    for i=1:p
        Omega(:,:,i)=Gamma'*sigRes(:,:,i)*Gamma;
        Sigma(:,:,i)=Gamma*Omega(:,:,i)*Gamma'+Gamma0*Omega0*Gamma0';
    end
    
    Yfit=zeros(n,r);
    for i=1:p
        if i>1
            Yfit(ind(ncum(i-1)+1:ncum(i)),:)=ones(ng(i),1)*mug(:,i)';
        else
            Yfit(ind(1:ncum(1)),:)=ones(ng(1),1)*mug(:,1)';
        end
    end

    fracN=ng/n;
    J=zeros(p*r+p*r*(r+1)/2);
    for i=1:p-1
        for j=1:p-1

            J((i-1)*r+1:i*r,(j-1)*r+1:j*r)=fracN(j)*fracN(i)/fracN(p)*inv(Sigma(:,:,p));
        end
        J((i-1)*r+1:i*r,(i-1)*r+1:i*r)=fracN(i)*inv(Sigma(:,:,i))+fracN(i)^2/fracN(p)*inv(Sigma(:,:,p));
    end
    for i=1:p
        J(r*(p-1)+1+(i-1)*r*(r+1)/2:r*(p-1)+i*r*(r+1)/2,r*(p-1)+1+(i-1)*r*(r+1)/2:r*(p-1)+i*r*(r+1)/2)=0.5*fracN(i)*Expan(r)'*kron(inv(Sigma(:,:,i)),inv(Sigma(:,:,i)))*Expan(r);
    end
    J(r+1:end,r+1:end)=J(1:(p-1)*r+p*r*(r+1)/2,1:(p-1)*r+p*r*(r+1)/2);
    J(1:r,:)=0;
    J(r+1:end,1:r)=0;
    for i=1:p
        J(1:r,1:r)=J(1:r,1:r)+fracN(i)*inv(Sigma(:,:,i));
    end
    for i=1:p-1
        J(1:r,i*r+1:(i+1)*r)=fracN(i)*(inv(Sigma(:,:,p))-inv(Sigma(:,:,i)));
        J(i*r+1:(i+1)*r,1:r)=J(1:r,i*r+1:(i+1)*r);
    end
    temp=inv(J);
    asyFm=reshape(sqrt(diag(temp(1:r*p,1:r*p))),r,p);
    
    H=zeros(p*r+p*r*(r+1)/2,r+u*(r+p-1-u)+p*u*(u+1)/2+(r-u)*(r-u+1)/2);
    for i=1:p-1
        H((i-1)*r+1:i*r,(i-1)*u+1:i*u)=Gamma;
        H((i-1)*r+1:i*r,(p-1)*u+1:u*(r+p-1-u))=kron(eta(:,i)',Gamma0);
    end
    for i=1:p 
        H(r*(p-1)+(i-1)*r*(r+1)/2+1:r*(p-1)+i*r*(r+1)/2,(p-1)*u+1:u*(r+p-1-u))=2*Contr(r)*(kron(Gamma*Omega(:,:,i),Gamma0)-kron(Gamma,Gamma0*Omega0));
        H(r*(p-1)+(i-1)*r*(r+1)/2+1:r*(p-1)+i*r*(r+1)/2,u*(r+p-1-u)+(i-1)*u*(u+1)/2+1:u*(r+p-1-u)+i*u*(u+1)/2)=Contr(r)*kron(Gamma,Gamma)*Expan(u);
        H(r*(p-1)+(i-1)*r*(r+1)/2+1:r*(p-1)+i*r*(r+1)/2,u*(r+p-1-u)+p*u*(u+1)/2+1:u*(r+p-1-u)+p*u*(u+1)/2+(r-u)*(r-u+1)/2)=Contr(r)*kron(Gamma0,Gamma0)*Expan(r-u);
    end
    H(r+1:end,r+1:end)=H(1:(p-1)*r+p*r*(r+1)/2,1:u*(r+p-1-u)+p*u*(u+1)/2+(r-u)*(r-u+1)/2);
    H(1:r,r+1:end)=0;
    H(r+1:end,1:r)=0;
    H(1:r,1:r)=eye(r);
    temp=H*inv(H'*J*H)*H';
    asyHenv=reshape(sqrt(diag(temp(1:r*p,1:r*p))),r,p);

    stat.mu=mu;
    stat.mug=mug;
    stat.Yfit=Yfit;
    stat.Gamma=Gamma;
    stat.Gamma0=Gamma0;
    stat.beta=beta;
    stat.Sigma=Sigma;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;
    stat.np=(r-u)+u*(r-u+p)+p*u*(u+1)/2+(r-u)*(r-u+1)/2;
    stat.l=-n*r/2*(1+log(2*pi))-n/2*(l+log(prod(eigtem(eigtem>0))));   
    stat.asyHenv=asyHenv;
    stat.ratio=asyFm./asyHenv;
    
end