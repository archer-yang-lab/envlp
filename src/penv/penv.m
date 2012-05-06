%% penv
% Fit the partial envelope model.

%% Usage
% stat = penv(X1,X2,Y,u,opts)
%
% Input
%
% X1: Predictors of main interst. An n by p1 matrix, n is the number of 
% observations, and p1 is the number of main predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
%
% X2: Covariates, or predictors not of main interest.  An n by p2 matrix,
% p2 is the number of covariates.  The covariates can be univariate or 
% multivariate, discrete or continuous.
%
% Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p1.
%
% u: Dimension of the partial envelope. An integer between 0 and r.
%
% opts: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * opts.maxIter: Maximum number of iterations.  Default value: 300.
% * opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * opts.verbose: Flag for print out output, logical 0 or 1. Default value:
% 0.
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
% * stat.covMatrix: The asymptotic covariance of (vec($$\beta_2$)', vec($$\beta_1$)')'.  An rp by
% rp matrix.  The covariance matrix returned are asymptotic.  For the
% actual standard errors, multiply by 1/n.
% * stat.asyPenv: Asymptotic standard error for elements in $$\beta_1$ under
% the partial envelope model.  An r by p1 matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1/sqrt(n).
% * stat.ratio: The asymptotic standard error ratio of the stanard multivariate 
% linear regression estimator over the partial envelope estimator, for each element 
% in $$\beta_1$.  An r by p1 matrix.
% * stat.np: The number of parameters in the envelope model.  A positive
% integer.
% * stat.n: The number of observations in the data.  A positive
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
% The following codes reconstruct the results of the paper and fiber
% example in Su and Cook (2012).
% 
% load T7-7.dat
% Y = T7_7(:,1:4);
% X = T7_7(:,5:7);
% X1 = X(:,3);
% X2 = X(:,1:2);
% alpha = 0.01;
% u = lrt_penv(X1,X2,Y,alpha)
% stat = penv(X1,X2,Y,u)
% stat.Omega
% eig(stat.Omega0)
% stat.ratio

function stat = penv(X1,X2,Y,u,opts)

% Verify and initialize the parameters

if (nargin < 4)
    error('Inputs: X1, X2, Y and u should be specified!');
elseif (nargin==4)
    opts = [];
end

[n,p1] = size(X1);
[n2,p2] = size(X2);
[n1,r] = size(Y);

if (n ~= n1||n2 ~= n1)
    error('The number of observations in X1, X2 and Y should be equal!');
end

if (p1 >= r)
    error('When the number of responses is less than the number of main predictors, the partial envelope model cannot be applied.');
end

u = floor(u);
if (u < 0 || u > r)
    error('u should be an integer between [0, r]!');
end

opts = make_opts(opts);

if isfield(opts,'init')
    [r2,u2] = size(opts.init);

    if (r ~= r2 || u ~= u2)
        error('The size of the initial value should be r by u!');
    end

    if (rank(opts.init) < u2)
        error('The initial value should be full rank!');
    end
end

%---preparation---


X1C = center(X1);
X2C = center(X2);
YC = center(Y);

SX1 = cov(X1,1);
SX2 = cov(X2,1);
SX12 = X1C'*X2C/n;

QX2 = eye(n)-X2C*inv(X2C'*X2C)*X2C';
R12 = QX2*X1C;
RY2 = QX2*YC;


% With different u, the model will be different.  When u = 0, X and Y are
% uncorrelated, so it should be fitted differently.  When u = r, the partial 
% envelope model reduces to the standard model, and it also should be fitted
% differently.


if u>0 && u<r


    %---Call env to compute most of the components in output---

    temp = env(R12,RY2,u,opts);
    beta1 = temp.beta;
    
    Gamma = temp.Gamma;
    Gamma0 = temp.Gamma0;
    eta = temp.eta;
    Omega = temp.Omega;
    Omega0 = temp.Omega0;
    Sigma1 = Gamma*Omega*Gamma';
    Sigma = temp.Sigma;
    
    
    stat.beta1 = beta1;
    stat.Gamma = Gamma;
    stat.eta = eta;
    stat.Omega = Omega;
    stat.Omega0 = Omega0;
    stat.Sigma = Sigma;
    stat.l = temp.l; 
    

    %---Compute the rest in output---
    beta2 = (YC-X1C*beta1')'*X2C*inv(X2C'*X2C);
    alpha = mean(Y)'-beta1*mean(X1)'-beta2*mean(X2)';
    
    stat.beta2 = beta2;    
    stat.alpha = alpha;
    stat.np = r+u*p1+r*p2+r*(r+1)/2;
    
    %---compute asymptotic variance and get the ratios---
    Sig1G2 = SX1-SX12*inv(SX2)*SX12';
    asyFm = kron(inv(Sig1G2),Sigma);
    asyFm = reshape(sqrt(diag(asyFm)),r,p1);
    
    temp=kron(eta*Sig1G2*eta',inv(Omega0))+kron(Omega,inv(Omega0))+kron(inv(Omega),Omega0)-2*kron(eye(u),eye(r-u));
    asyPenv=kron(inv(Sig1G2),Sigma1)+kron(eta',Gamma0)*inv(temp)*kron(eta,Gamma0');
    asyPenv=reshape(sqrt(diag(asyPenv)),r,p1);
    
    p=p1+p2;
    invSigma=inv(Sigma);
    J=zeros(r*p+r*(r+1)/2);
    J(1:p2*r,1:p2*r)=kron(SX2,invSigma);
    J(1+p2*r:p*r,1+p2*r:p*r)=kron(SX1,invSigma);
    J(1+p2*r:p*r,1:p2*r)=kron(SX12,invSigma);
    J(1:p2*r,1+p2*r:p*r)=J(1+p2*r:p*r,1:p2*r)';
    J(r*p+1:end,r*p+1:end)=Expan(r)'*kron(invSigma,invSigma)*Expan(r)/2;
    H=zeros(r*p+r*(r+1)/2,r*p2+u*p1+r*(r+1)/2);
    H(1:r*p2,1:r*p2)=eye(r*p2);
    H(r*p2+1:r*p,r*p2+1:r*p2+u*p1)=kron(eye(p1),Gamma);
    H(r*p2+1:r*p,r*p2+u*p1+1:r*p2+u*p1+u*(r-u))=kron(eta',Gamma0);
    H(r*p+1:end,r*p2+u*p1+1:r*p2+u*p1+u*(r-u)) = 2*Contr(r)*(kron(Gamma*Omega,Gamma0)-kron(Gamma,Gamma0*Omega0));
    H(r*p+1:end,r*p2+u*p1+u*(r-u)+1:r*p2+u*p1+u*(r-u)+u*(u+1)/2) = Contr(r)*kron(Gamma,Gamma)*Expan(u);
    H(r*p+1:end,r*p2+u*p1+u*(r-u)+u*(u+1)/2+1:end) = Contr(r)*kron(Gamma0,Gamma0)*Expan(r-u);
    covMatrix = H*inv(H'*J*H)*H';
    covMatrix = covMatrix(1:r*p,1:r*p);
    
    stat.covMatrix = covMatrix;
    stat.asyPenv = asyPenv;
    stat.ratio = asyFm./asyPenv;
    stat.n = n;
    
elseif u==0
    
    temp = fit_OLS(X2,Y);
    beta2 = temp.betaOLS;
    Sigma = temp.SigmaOLS;
    eigtem = eig(Sigma);
    
    stat.beta1 = zeros(r,p1);
    stat.beta2 = beta2;
    stat.Sigma = Sigma;
    stat.eta = [];
    stat.Gamma = [];
    stat.Gamma0 = eye(r);
    stat.Omega = [];
    stat.Omega0 = Sigma;
    stat.alpha = mean(Y)'-beta2*mean(X2)';
    stat.l = -n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.covMatrix = [];
    stat.asyPenv = [];
    stat.ratio = ones(r,p1);
    stat.np = r+u*p1+r*p2+r*(r+1)/2;
    stat.n = n;    

elseif u==r
    
    X = [X1 X2];
    temp = fit_OLS(X,Y);
    beta = temp.betaOLS;
    Sigma = temp.SigmaOLS;
    eigtem = eig(Sigma);
    sigX = cov(X,1);
    tempasy = kron(inv(sigX),Sigma);
    covMatrix = tempasy(1:r*p1,1:r*p1);
    asyFm = reshape(sqrt(diag(covMatrix)),r,p1);
    
    
    stat.beta1 = beta(:,1:p1);
    stat.beta2 = beta(:,p1+1:end);
    stat.Sigma = Sigma;
    stat.eta = beta(:,1:p1);
    stat.Gamma = eye(r);
    stat.Gamma0 = [];
    stat.Omega = Sigma;
    stat.Omega0 = [];
    stat.alpha = mean(Y)'-beta*mean(X)';
    stat.l = -n*r/2*(1+log(2*pi))-n/2*log(prod(eigtem(eigtem>0)));
    stat.covMatrix = covMatrix;
    stat.asyPenv = asyFm;
    stat.ratio = ones(r,p1);
    stat.np = r+u*p1+r*p2+r*(r+1)/2;
    stat.n = n;
    
end
    