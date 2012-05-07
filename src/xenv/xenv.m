%% xenv
% Fit the envelope model for the reduction on X.

%% Syntax
% ModelOutput=xenv(X,Y,u)
% ModelOutput=xenv(X,Y,u,Opts)
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
% Opts: A list containing the optional input parameter, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out output, logical 0 or 1. Default value:
% 0. 
%
% Output
% 
% ModelOutput: A list that contains the maximum likelihood estimators and some
% statistics.
% 
% * ModelOutput.beta: The envelope estimator of the regression coefficients $$\beta$. 
% An p by r matrix.
% * ModelOutput.SigX: The envelope estimator of the covariance matrix of X, $$\Sigma_X$.  A p by
% p matrix.
% * ModelOutput.Gamma: The orthogonal basis of the envelope subspace. An p by u
% semi-orthogonal matrix.
% * ModelOutput.Gamma0: The orthogonal basis of the complement of the envelope
% subspace.  An p by p-u semi-orthogonal matrix.
% * ModelOutput.eta: The coordinates of $$\beta$ with respect to Gamma. An u by r
% matrix.
% * ModelOutput.Omega: The coordinates of $$\Sigma_X$ with respect to Gamma. An u by u
% matrix.
% * ModelOutput.Omega0: The coordinates of $$\Sigma_X$ with respect to Gamma0. An p-u by p-u
% matrix.
% * ModelOutput.mu: The estimated intercept.  An r by 1 vector.
% * ModelOutput.sigYcX: The estimated conditional covariance matrix of Y given X.
% An r by r matrix.
% * ModelOutput.l: The maximized log likelihood function.  A real number.
% * ModelOutput.covMatrix: The asymptotic covariance of vec($$\beta$).  An pr by
% pr matrix.  The covariance matrix returned are asymptotic.  For the
% actual standard errors, multiply by 1/n.
% * ModelOutput.asyXenv: Asymptotic standard error for elements in $$\beta$ under
% the envelope model.  An r by p matrix.  The standard errors returned are
% asymptotic, for actual standard errors, multiply by 1/sqrt(n).
% * ModelOutput.ratio: The asymptotic standard error ratio of the standard multivariate 
% linear regression estimator over the envelope estimator, for each element 
% in $$\beta$.  An p by r matrix.
% * ModelOutput.np: The number of parameters in the envelope model.  A positive
% integer.
% * ModelOutput.n: The number of observations in the data.  A positive
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
% ModelOutput=xenv(X,Y,0);
% 
% p=size(X,2);
% ModelOutput=xenv(X,Y,p);
% 
% % When u=p, the envelope model reduces to the ordinary least squares
% % regression
% 
% temp=fit_OLS(X,Y);
% temp.SigmaOLS
% ModelOutput.sigYcX
% temp.betaOLS'
% ModelOutput.beta
% 
% ModelOutput=xenv(X,Y,5);
% 
% %  To compare with the results obtained by Partial Least Squares, use the
% %  command 
% % [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,5);
% % ModelOutput.beta
% % ModelOutput.mu
% % BETA

function ModelOutput=xenv(X,Y,u,Opts)

if (nargin < 3)
    error('Inputs: X, Y and u should be specified!');
elseif (nargin==3)
    Opts=[];
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

Opts=make_opts(Opts);

if isfield(Opts,'init')
    [r2,u2]=size(Opts.init);

    if (r ~= r2 || u ~= u2)
        error('The size of the initial value should be r by u!');
    end

    if (rank(Opts.init) < u2)
        error('The initial value should be full rank!');
    end
end

%---preparation---
DataParameter=make_parameter(X,Y,'xenv');

n=DataParameter.n;
p=DataParameter.p;
r=DataParameter.r;
mX=DataParameter.mX;
mY=DataParameter.mY;
sigX=DataParameter.sigX; % Standard estimator of \Sigma_X
sigXY=DataParameter.sigXY;
sigY=DataParameter.sigY;
invSigX=DataParameter.invSigX;

sigYcX=sigY-sigXY'*invSigX*sigXY;

% With different u, the model will be different.  When u=0, X and Y are
% uncorrelated, so it should be fitted differently.  When u=r, the envelope
% model reduces to the standard model, and it also should be fitted
% differently.


if u>0 && u<p


    %---Compute \Gamma using sg_min---
    tempParameter=make_parameter(Y,X,'env');
    tempF = make_F(@F4env,tempParameter);
    
    
    maxIter=Opts.maxIter;
	ftol=Opts.ftol;
	gradtol=Opts.gradtol;
	if (Opts.verbose==0) 
        verbose='quiet';
    else
        verbose='verbose';
    end
    if ~isfield(Opts,'init') 
        init=get_Init(tempF,Y,X,u,tempParameter);
    else
        init=Opts.init;
    end
    
    
    F = make_F(@F4xenv,DataParameter);
    dF = make_dF(@dF4xenv,DataParameter);

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
    

    ModelOutput.beta=beta;
    ModelOutput.SigX=SigX;
    ModelOutput.Gamma=Gamma;
    ModelOutput.Gamma0=Gamma0;
    ModelOutput.eta=eta;
    ModelOutput.Omega=Omega;
    ModelOutput.Omega0=Omega0;
    ModelOutput.mu=mu;
    ModelOutput.sigYcX=sigYcX;
    ModelOutput.l=-n*(p+r)/2*log(2*pi)-n/2*(log(prod(eigtem(eigtem>0)))+log(prod(eigtem2(eigtem2>0))))-1/2*trace(X*invSigX*X')-1/2*trace((Y-ones(n,1)*mu'-X*beta)*inv(sigYcX)*(Y-ones(n,1)*mu'-X*beta)');
    ModelOutput.covMatrix=covMatrix;
    ModelOutput.asyXenv=asyEnv;
    ModelOutput.ratio=asyFm./asyEnv;
    ModelOutput.np=r+u*r+p*(p+1)/2+r*(r+1)/2;
    ModelOutput.n=n;    
    
elseif u==0
    
    mu=mY;
    eigtem=eig(sigX);
    eigtem3=eig(sigY);    
    ModelOutput.beta=zeros(p,r);
    ModelOutput.SigX=sigX;
    ModelOutput.Gamma=[];
    ModelOutput.Gamma0=eye(p);
    ModelOutput.eta=[];
    ModelOutput.Omega=[];
    ModelOutput.Omega0=sigX;
    ModelOutput.mu=mu;
    ModelOutput.sigYcX=sigY;
    ModelOutput.l=-n*(p+r)/2*log(2*pi)-n/2*(log(prod(eigtem(eigtem>0)))+log(prod(eigtem3(eigtem3>0))))-1/2*trace(X*invSigX*X')-1/2*trace((Y-ones(n,1)*mu')*inv(sigY)*(Y-ones(n,1)*mu')');
    ModelOutput.covMatrix=[];
    ModelOutput.asyXenv=[];
    ModelOutput.ratio=ones(p,r);
    ModelOutput.np=r+p*(p+1)/2+r*(r+1)/2;
    ModelOutput.n=n;    

elseif u==p
    
    temp=fit_OLS(X,Y);
 
    covMatrix=kron(sigYcX,invSigX);
    asyFm=reshape(sqrt(diag(covMatrix)),p,r);
    beta=temp.betaOLS';
    mu=mY-beta'*mX;
    eigtem=eig(sigX);
    eigtem2=eig(sigYcX);
    
    ModelOutput.beta=beta;
    ModelOutput.SigX=sigX;
    ModelOutput.Gamma=eye(p);
    ModelOutput.Gamma0=[];
    ModelOutput.eta=beta;
    ModelOutput.Omega=sigX;
    ModelOutput.Omega0=[];
    ModelOutput.mu=mu;
    ModelOutput.sigYcX=sigYcX;
    ModelOutput.l=-n*(p+r)/2*log(2*pi)-n/2*(log(prod(eigtem(eigtem>0)))+log(prod(eigtem2(eigtem2>0))))-1/2*trace(X*invSigX*X')-1/2*trace((Y-ones(n,1)*mu'-X*beta)*inv(sigYcX)*(Y-ones(n,1)*mu'-X*beta)');
    ModelOutput.covMatrix=covMatrix;
    ModelOutput.asyXenv=asyFm;
    ModelOutput.ratio=ones(p,r);
    ModelOutput.np=r+(p+r)*(p+r+1)/2;
    ModelOutput.n=n;    
    
end
    
    
    