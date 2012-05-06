%% xenv
% Fit the envelope model for the reduction on X.

%% Syntax
% stat=xenv(X,Y,u)
% stat=xenv(X,Y,u,opts)
%
% Input
%
% X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
%
% Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p.
%
% u: Dimension of the envelope. An integer between 0 and p.
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
% * stat.beta: The envelope estimator of the regression coefficients $$\beta$. 
% An p by r matrix.
% * stat.SigX: The envelope estimator of the covariance matrix of X, $$\Sigma_X$.  A p by
% p matrix.
% * stat.Gamma: The orthogonal basis of the envelope subspace. An p by u
% semi-orthogonal matrix.
% * stat.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An p by p-u semi-orthogonal matrix.
% * stat.eta: The coordinates of $$\beta$ with respect to Gamma. An u by r
% matrix.
% * stat.Omega: The coordinates of $$\Sigma_X$ with respect to Gamma. An u by u
% matrix.
% * stat.Omega0: The coordinates of $$\Sigma_X$ with respect to Gamma0. An p-u by p-u
% matrix.
% * stat.mu: The estimated intercept.  An r by 1 vector.
% * stat.sigYcX: The estimated conditional covariance matrix of Y given X.
% An r by r matrix.
% * stat.l: The maximized log likelihood function.  A real number.
% * stat.covMatrix: The asymptotic covariance of vec($$\beta$).  An pr by
% pr matrix.  The covariance matrix returned are asymptotic.  For the
% actual standard errors, multiply by 1/n.
% * stat.asyXenv: Asymptotic standard error for elements in $$\beta$ under
% the envelope model.  An r by p matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1/sqrt(n).
% * stat.ratio: The asymptotic standard error ratio of the standard multivariate 
% linear regression estimator over the envelope estimator, for each element 
% in $$\beta$.  An p by r matrix.
% * stat.np: The number of parameters in the envelope model.  A positive
% integer.
% * stat.n: The number of observations in the data.  A positive
% integer.

%% Description
% This function fits the envelope model to the responses and predictors,
% using the maximum likehood estimation.  When the dimension of the
% envelope is between 1 and r-1, we implemented the algorithm in Cook et
% al. (2012).  When the dimension is r, then the envelope model degenerates
% to the standard multivariate linear regression.  When the dimension is 0,
% it means that X and Y are uncorrelated, and the fitting is different.

%% References
% 
% * The codes is implemented based on the algorithm in Section 4.5.1 of Cook 
% et al (2012).
% * The Grassmann manifold optimization step calls the package sg_min 2.4.1
% by Ross Lippert (http://web.mit.edu/~ripper/www.sgmin.html).

%% Example
%
% load wheatprotein.txt
% X=wheatprotein(:,1:6);
% Y=wheatprotein(:,7);
% stat=xenv(X,Y,0);
% 
% p=size(X,2);
% stat=xenv(X,Y,p);
% 
% % When u=p, the envelope model reduces to the ordinary least squares
% % regression
% 
% temp=fit_OLS(X,Y);
% temp.SigmaOLS
% stat.sigYcX
% temp.betaOLS'
% stat.beta
% 
% stat=xenv(X,Y,5);
% 
% %  To compare with the results obtained by Partial Least Squares, use the
% %  command 
% % [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,5);
% % stat.beta
% % stat.mu
% % BETA

function stat=xenv(X,Y,u,opts)

if (nargin < 3)
    error('Inputs: X, Y and u should be specified!');
elseif (nargin==3)
    opts=[];
end

[n,p]=size(X);
[n1,r]=size(Y);

if (n ~= n1)
    error('The number of observations in X and Y should be equal!');
end

if (r >= p)
    error('When the number of predictors is less than the number of responses, the envelope model for reduction on X cannot be applied.');
end

u = floor(u);
if (u < 0 || u > p)
    error('u should be an integer between [0, p]!');
end

opts=make_opts(opts);

if isfield(opts,'init')
    [r2,u2]=size(opts.init);

    if (r ~= r2 || u ~= u2)
        error('The size of the initial value should be r by u!');
    end

    if (rank(opts.init) < u2)
        error('The initial value should be full rank!');
    end
end

%---preparation---
dataParameter=make_parameter(X,Y,'xenv');

n=dataParameter.n;
p=dataParameter.p;
r=dataParameter.r;
mX=dataParameter.mX;
mY=dataParameter.mY;
sigX=dataParameter.sigX; % Standard estimator of \Sigma_X
sigXY=dataParameter.sigXY;
sigY=dataParameter.sigY;
invSigX=dataParameter.invSigX;

sigYcX=sigY-sigXY'*invSigX*sigXY;

% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the envelope
% model reduces to the standard model, and it also should be fitted
% differently.


if u>0 && u<p


    %---Compute \Gamma using sg_min---
    tempParameter=make_parameter(Y,X,'env');
    tempF = make_F(@F4env,tempParameter);
    
    
    maxIter=opts.maxIter;
	ftol=opts.ftol;
	gradtol=opts.gradtol;
	if (opts.verbose==0) 
        verbose='quiet';
    else
        verbose='verbose';
    end
    if ~isfield(opts,'init') 
        init=get_Init(tempF,Y,X,u,tempParameter);
    else
        init=opts.init;
    end
    
    
    F = make_F(@F4xenv,dataParameter);
    dF = make_dF(@dF4xenv,dataParameter);

    [l Gamma]=sg_min(F,dF,init,maxIter,'prcg',verbose,ftol,gradtol);



    %---Compute the rest of the parameters based on \Gamma---
    Gamma0=grams(nulbasis(Gamma'));
    eta=Gamma'*sigXY;
    Omega=Gamma'*sigX*Gamma;
    Omega0=Gamma0'*sigX*Gamma0;
    SigX=Gamma*Omega*Gamma'+Gamma0*Omega0*Gamma0'; % Envelope estimator of \Sigma_X
    beta=Gamma*inv(Omega)*eta;
    mu=mY-beta'*mX;

    eigtem=eig(SigX);
    eigtem2=eig(sigYcX);

    %---compute asymptotic variance and get the ratios---
    asyFm=kron(sigYcX,inv(SigX));
    asyFm=reshape(sqrt(diag(asyFm)),p,r);
    temp=kron(eta*inv(sigYcX)*eta',Omega0)+kron(Omega,inv(Omega0))+kron(inv(Omega),Omega0)-2*kron(eye(u),eye(p-u));
    covMatrix=kron(sigYcX,Gamma*inv(Omega)*Gamma')+kron(eta',Gamma0)*inv(temp)*kron(eta,Gamma0');
    asyEnv=reshape(sqrt(diag(covMatrix)),p,r);
    

    stat.beta=beta;
    stat.SigX=SigX;
    stat.Gamma=Gamma;
    stat.Gamma0=Gamma0;
    stat.eta=eta;
    stat.Omega=Omega;
    stat.Omega0=Omega0;
    stat.mu=mu;
    stat.sigYcX=sigYcX;
    stat.l=-n*(p+r)/2*log(2*pi)-n/2*(log(prod(eigtem(eigtem>0)))+log(prod(eigtem2(eigtem2>0))))-1/2*trace(X*invSigX*X')-1/2*trace((Y-ones(n,1)*mu'-X*beta)*inv(sigYcX)*(Y-ones(n,1)*mu'-X*beta)');
    stat.covMatrix=covMatrix;
    stat.asyXenv=asyEnv;
    stat.ratio=asyFm./asyEnv;
    stat.np=r+u*r+p*(p+1)/2+r*(r+1)/2;
    stat.n=n;    
    
elseif u==0
    
    mu=mY;
    eigtem=eig(sigX);
    eigtem3=eig(sigY);    
    stat.beta=zeros(p,r);
    stat.SigX=sigX;
    stat.Gamma=[];
    stat.Gamma0=eye(p);
    stat.eta=[];
    stat.Omega=[];
    stat.Omega0=sigX;
    stat.mu=mu;
    stat.sigYcX=sigY;
    stat.l=-n*(p+r)/2*log(2*pi)-n/2*(log(prod(eigtem(eigtem>0)))+log(prod(eigtem3(eigtem3>0))))-1/2*trace(X*invSigX*X')-1/2*trace((Y-ones(n,1)*mu')*inv(sigY)*(Y-ones(n,1)*mu')');
    stat.covMatrix=[];
    stat.asyXenv=[];
    stat.ratio=ones(p,r);
    stat.np=r+p*(p+1)/2+r*(r+1)/2;
    stat.n=n;    

elseif u==p
    
    temp=fit_OLS(X,Y);
 
    covMatrix=kron(sigYcX,invSigX);
    asyFm=reshape(sqrt(diag(covMatrix)),p,r);
    beta=temp.betaOLS';
    mu=mY-beta'*mX;
    eigtem=eig(sigX);
    eigtem2=eig(sigYcX);
    
    stat.beta=beta;
    stat.SigX=sigX;
    stat.Gamma=eye(p);
    stat.Gamma0=[];
    stat.eta=beta;
    stat.Omega=sigX;
    stat.Omega0=[];
    stat.mu=mu;
    stat.sigYcX=sigYcX;
    stat.l=-n*(p+r)/2*log(2*pi)-n/2*(log(prod(eigtem(eigtem>0)))+log(prod(eigtem2(eigtem2>0))))-1/2*trace(X*invSigX*X')-1/2*trace((Y-ones(n,1)*mu'-X*beta)*inv(sigYcX)*(Y-ones(n,1)*mu'-X*beta)');
    stat.covMatrix=covMatrix;
    stat.asyXenv=asyFm;
    stat.ratio=ones(p,r);
    stat.np=r+(p+r)*(p+r+1)/2;
    stat.n=n;    
    
end
    
    
    