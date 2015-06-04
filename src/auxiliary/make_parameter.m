%% make_parameter
% Compute summary statistics from the data.

%% Syntax
%         DataParameter = make_parameter(X, Y, method)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors. The
% predictors can be univariate or multivariate, discrete or continuous.
% 
% *Y*: Multivariate responses. An n by r matrix, r is the number of
% responses and n is number of observations. The responses must be 
% continuous variables, and r should be strictly greater than p.
% 
% *method*: A string of characters indicating which member of the envelope
% family to be used, the choices can be 'env', 'ienv', 'henv', 'senv', 
% 'sxenv', 'xenv' or 'xenvpls'.
%
%% Output
% 
% *DataParameter*: A list that contains summary statistics computed from the
% data.  The output list can vary from method to method.
% 
% * DataParameter.n: The number of observations in the data.  A positive
% integer.  
% * DataParameter.ng: A p by 1 vector containing the number of observations
% in each group.  p is the number of groups.  Only for 'henv'.
% * DataParameter.ncum: A p by 1 vector containing the total number of
% observations till this group.  Only for 'henv'.
% * DataParameter.ind: An n by 1 vector indicating the sequence of the
% observations after sorted by groups.
% * DataParameter.p: The number of predictors or number of groups for
% 'henv'.  A positive integer.
% * DataParameter.r: The number of responses.  A positive integer.
% * DataParameter.XC: Centered predictors.  An n by p matrix with the ith
% row being the ith observation of X subtracted by the mean of X.  Only for
% 'env' and 'ienv'.
% * DataParameter.YC: Centered responses.  An n by r matrix with the ith
% row being the ith observation of Y subtracted by the mean of Y.  Only for
% 'env' and 'ienv'.
% * DataParameter.mX: The mean of predictors.  A p by 1 vector.  For all
% method except 'henv'.
% * DataParameter.mY: The mean of responses.  An r by 1 vector.
% * DataParameter.mYg: An r by p matrix with the ith column being the
% sample mean of the ith group.
% * DataParameter.sigX: The sample covariance matrix of X.  A p by p
% matrix.
% * DataParameter.sigY: The sample covariance matrix of Y.  An r by r
% matrix.
% * DataParameter.sigRes: For 'env', 'senv', 'ienv': The sample covariance
% matrix of the residuals from the ordinary least squares regression of Y
% on X.  An r by r matrix. For 'henv', an r by r by p three dimensional
% matrix with the ith depth is the ith sample covariance matrix for the ith
% group.
% * DataParameter.sigFit: The sample covariance matrix of the fitted value 
% from the ordinary least squares regression of Y on X.  An r by r matrix.
% Only for method 'ienv'.
% * DataParameter.betaOLS: The regression coefficients from the ordinary
% least squares regression of Y on X.  An r by p matrix.  For all methods
% except 'henv'.
% * DataParameter.invsigY: The inverse of the sample covariance matrix of Y.  An r by r
% matrix. For all methods except 'ienv', 'xenv' and 'xenvpls'.
% * DataParameter.invsigRes: The inverse of the sample covariance matrix of 
% the residuals form the ordinary least squares regression of Y on X.  An r by r
% matrix. Only for method 'ienv'.

%% Description
% This function computes statistics that will be used frequently in the
% estimation for each method.


function DataParameter = make_parameter(X, Y, method)

if (strcmp(method, 'env'))
    
    
    [n, p] = size(X);
    r = size(Y, 2);

    XC = center(X);
    YC = center(Y);
    ModelOutput = fit_OLS(X, Y);
    sigY = cov(Y, 1);
    invsigY = inv(sigY);
    eigtemY = eig(sigY);
    logDetSigY = sum(log(eigtemY(eigtemY > 0)));

    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.XC = XC;
    DataParameter.YC = YC;
    DataParameter.mX = mean(X)';
    DataParameter.mY = mean(Y)';
    DataParameter.sigX = cov(X, 1);
    DataParameter.sigY = sigY;
    DataParameter.sigRes = ModelOutput.SigmaOLS;
    DataParameter.betaOLS = ModelOutput.betaOLS;
    DataParameter.logDetSigY = logDetSigY;
    DataParameter.invsigY = invsigY;
    
elseif (strcmp(method, 'senv'))
    
    
    [n, p] = size(X);
    r = size(Y, 2);

    ModelOutput = fit_OLS(X, Y);
    sigY = cov(Y, 1);
    invsigY = inv(sigY);
    eigtemY = eig(sigY);
    logDetSigY = sum(log(eigtemY(eigtemY > 0)));
    
    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.mX = mean(X)';
    DataParameter.mY = mean(Y)';
    DataParameter.sigX = cov(X, 1);
    DataParameter.sigY = sigY;
    DataParameter.sigRes = ModelOutput.SigmaOLS;
    DataParameter.betaOLS = ModelOutput.betaOLS;
    DataParameter.logDetSigY = logDetSigY;
    DataParameter.invsigY = invsigY;
    
    
elseif (strcmp(method, 'ienv'))
    
    
    [n, p] = size(X);
    r = size(Y, 2);

    XC = center(X);
    YC = center(Y);
    ModelOutput = fit_OLS(X, Y);
    sigY = cov(Y, 1);
    sigRes = ModelOutput.SigmaOLS;
    sigFit = sigY - sigRes;
    eigtem = eig(sigRes);
    logDetSigRes = sum(log(eigtem(eigtem > 0)));
    invsigRes = inv(sigRes);

    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.XC = XC;
    DataParameter.YC = YC;
    DataParameter.mX = mean(X)';
    DataParameter.mY = mean(Y)';
    DataParameter.sigX = cov(X, 1);
    DataParameter.sigY = sigY;
    DataParameter.sigRes = sigRes;
    DataParameter.sigFit = sigFit;
    DataParameter.betaOLS = ModelOutput.betaOLS;
    DataParameter.logDetSigRes = logDetSigRes;
    DataParameter.invsigRes = invsigRes;
    
    
elseif (strcmp(method, 'henv'))
    
    
    [n, r] = size(Y);

    [Xs, ind] = sortrows(X);
    p = 1;
    temp = Xs(1, :);
    ng = 0;
    for i = 1 : n
        if prod(double(Xs(i, :) == temp)) == 0            
            temp = Xs(i, :);
            ng(p) = i - 1 - sum(ng(1 : p - 1));
            ncum(p) = i - 1;
            p = p + 1;
        end
    end           
    ncum(p) = n;
    ng(p) = n - sum(ng(1 : p - 1));

    sigRes = zeros(r, r, p);
    mYg = zeros(r, p);

    Ys = Y(ind, :);

    for i = 1 : p
        if i > 1
            sigRes(:, :, i) = cov(Ys(ncum(i - 1) + 1 : ncum(i), :), 1);
            mYg(:, i) = mean(Ys(ncum(i - 1) + 1 : ncum(i), :))';
        else
            sigRes(:, :, i) = cov(Ys(1 : ncum(i), :), 1);
            mYg(:, i) = mean(Ys(1 : ncum(i), :))';
        end
    end

    sigY = cov(Y, 1);
    eigtemY = eig(sigY);
    logDetSigY = sum(log(eigtemY(eigtemY > 0)));
    invsigY = inv(sigY);

    DataParameter.n = n;
    DataParameter.ng = ng;
    DataParameter.ncum = ncum;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.ind = ind;
    DataParameter.mY = mean(Y)';
    DataParameter.mYg = mYg;
    DataParameter.sigY = sigY;
    DataParameter.sigRes = sigRes;
    DataParameter.logDetSigY = logDetSigY;
    DataParameter.invsigY = invsigY;
    
elseif (strcmp(method, 'sxenv'))

    [n, p] = size(X);
    r = size(Y, 2);
    XC = center(X);
    YC = center(Y);
    sigX = cov(X, 1);
    sigXY = XC' * YC / n;
    sigY = cov(Y, 1);
    
    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.XC = XC;
    DataParameter.YC = YC;
    DataParameter.mX = mean(X)';
    DataParameter.mY = mean(Y)';
    DataParameter.sigX = sigX;
    DataParameter.sigY = sigY;
    DataParameter.sigXY = sigXY;
    DataParameter.sigXcY = sigX - sigXY / sigY * sigXY';
    DataParameter.invSigX = inv(sigX);

elseif (strcmp(method, 'xenv'))
    
    [n, p] = size(X);
    r = size(Y, 2);
    XC = center(X);
    YC = center(Y);
    sigX = cov(X, 1);
    sigXY = XC' * YC / n;
    sigY = cov(Y, 1);
    eigtemY = eig(sigY);
    logDetSigY = sum(log(eigtemY(eigtemY > 0)));
    eigtemX = eig(sigX);
    logDetSigX = sum(log(eigtemX(eigtemX > 0)));
    
    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.mX = mean(X)';
    DataParameter.mY = mean(Y)';
    DataParameter.sigX = sigX;
    DataParameter.sigY = sigY;
    DataParameter.sigXY = sigXY;
    DataParameter.sigXcY = sigX - sigXY / sigY * sigXY';
    DataParameter.invSigX = inv(sigX);
    DataParameter.logDetSigY = logDetSigY;
    DataParameter.logDetSigX = logDetSigX;
    
elseif (strcmp(method, 'xenvpls'))
    
    [n, p] = size(X);
    r = size(Y, 2);
    XC = center(X);
    YC = center(Y);
    sigX = cov(X, 1);
    sigXY = XC' * YC / n;
    sigY = cov(Y, 1);
    
    DataParameter.n = n;
    DataParameter.p = p;
    DataParameter.r = r;
    DataParameter.mX = mean(X)';
    DataParameter.mY = mean(Y)';
    DataParameter.sigX = sigX;
    DataParameter.sigY = sigY;
    DataParameter.sigXY = sigXY;

end


