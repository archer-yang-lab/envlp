%% make_parameter
% Compute summary statistics from the data.

%% Usage
% dataParameter=make_parameter(X,Y,method)
%
% Input
%
% * X: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% * Y: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p.
% * method: A string of characters indicating which member of the envelope
% family to be used, the choices can be 'env', 'ienv', 'henv' or 'senv'.  
%
% Output
% 
% dataParameter: A list that contains summary statistics computed from the
% data.  The output list can vary from method to method.
% 
% * dataParameter.n: The number of observations in the data.  A positive
% integer.  
% * dataParameter.ng: A p by 1 vector containing the number of observations
% in each group.  p is the number of groups.  Only for 'henv'.
% * dataParameter.ncum: A p by 1 vector containing the total number of
% observations till this group.  Only for 'henv'.
% * dataParameter.ind: An n by 1 vector indicating the sequence of the
% observations after sorted by groups.
% * dataParameter.p: The number of predictors or number of groups for
% 'henv'.  A positive integer.
% * dataParameter.r: The number of responses.  A positive integer.
% * dataParameter.XC: Centered predictors.  An n by p matrix with the ith
% row being the ith observation of X subtracted by the mean of X.  Only for
% 'env' and 'ienv'.
% * dataParameter.YC: Centered responses.  An n by r matrix with the ith
% row being the ith observation of Y subtracted by the mean of Y.  Only for
% 'env' and 'ienv'.
% * dataParameter.mX: The mean of predictors.  A p by 1 vector.  For all
% method except 'henv'.
% * dataParameter.mY: The mean of responses.  An r by 1 vector.
% * dataParameter.mYg: An r by p matrix with the ith column being the
% sample mean of the ith group.
% * dataParameter.sigX: The sample covariance matrix of X.  A p by p
% matrix.
% * dataParameter.sigY: The sample covariance matrix of Y.  An r by r
% matrix.
% * dataParameter.sigRes: For 'env', 'senv', 'ienv': The sample covariance
% matrix of the residuals from the ordinary least squares regression of Y
% on X.  An r by r matrix. For 'henv', an r by r by p three dimensional
% matrix with the ith depth is the ith sample covariance matrix for the ith
% group.
% * dataParameter.sigFit: The sample covariance matrix of the fitted value 
% from the ordinary least squares regression of Y on X.  An r by r matrix.
% Only for method 'ienv'.
% * dataParameter.betaOLS: The regression coefficients from the ordinary
% least squares regression of Y on X.  An r by p matrix.  For all method
% except 'henv'.

%% Description
% This function computes statistics that will be used frequently in the
% estimation for each method.


function dataParameter=make_parameter(X,Y,method)

if (strcmp(method,'env'))
    
    
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
   
    
elseif (strcmp(method,'senv'))
    
    
    [n p]=size(X);
    r=size(Y,2);

    XC=center(X);
    YC=center(Y);
    [betaOLS sigRes]=fit_OLS(X,Y);

    dataParameter.n=n;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.mX=mean(X)';
    dataParameter.mY=mean(Y)';
    dataParameter.sigX=cov(X,1);
    dataParameter.sigY=cov(Y,1);
    dataParameter.sigRes=sigRes;
    dataParameter.betaOLS=betaOLS;
    
    
elseif (strcmp(method,'ienv'))
    
    
    [n p]=size(X);
    r=size(Y,2);

    XC=center(X);
    YC=center(Y);
    [betaOLS sigRes]=fit_OLS(X,Y);
    sigY=cov(Y,1);
    sigFit=sigY-sigRes;

    dataParameter.n=n;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.XC=XC;
    dataParameter.YC=YC;
    dataParameter.mX=mean(X)';
    dataParameter.mY=mean(Y)';
    dataParameter.sigX=cov(X,1);
    dataParameter.sigY=sigY;
    dataParameter.sigRes=sigRes;
    dataParameter.sigFit=sigFit;
    dataParameter.betaOLS=betaOLS;
    
    
elseif (strcmp(method,'henv'))
    
    
    [n r]=size(Y);

    [Xs ind]=sortrows(X);
    p=1;
    temp=Xs(1,:);
    ng=0;
    for i=1:n
        if prod(double(Xs(i,:)==temp))==0            
            temp=Xs(i,:);
            ng(p)=i-1-sum(ng(1:p-1));
            ncum(p)=i-1;
            p=p+1;
        end
    end           
    ncum(p)=n;
    ng(p)=n-sum(ng(1:p-1));

    sigRes=zeros(r,r,p);
    mYg=zeros(r,p);

    Ys=Y(ind,:);

    for i=1:p
        if i>1
            sigRes(:,:,i)=cov(Ys(ncum(i-1)+1:ncum(i),:),1);
            mYg(:,i)=mean(Ys(ncum(i-1)+1:ncum(i),:))';
        else
            sigRes(:,:,i)=cov(Ys(1:ncum(i),:),1);
            mYg(:,i)=mean(Ys(1:ncum(i),:))';
        end
    end



    dataParameter.n=n;
    dataParameter.ng=ng;
    dataParameter.ncum=ncum;
    dataParameter.p=p;
    dataParameter.r=r;
    dataParameter.ind=ind;
    dataParameter.mY=mean(Y)';
    dataParameter.mYg=mYg;
    dataParameter.sigY=cov(Y,1);
    dataParameter.sigRes=sigRes;
    
    

end


