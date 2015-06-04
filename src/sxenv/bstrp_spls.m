%% bstrp_spls
% Compute bootstrap standard error for the scaled SIMPLS algorithm. 

%% Syntax
%         bootse = bstrp_spls(X, Y, u, B)
%         bootse = bstrp_spls(X, Y, u, B, Opts)
%
%% Input
%
% *X*: Predictors. An n by p matrix, p is the number of predictors and n is
% number of observations.  The predictors must be continuous variables.
% 
% *Y*: Responses. An n by r matrix, r is the number of
% responses. The response can be univariate or multivariate and must be
% continuous variable.
% 
% *u*: Dimension of the envelope subspace.  A positive integer between 0 and
% r.
% 
% *B*: Number of bootstrap samples.  A positive integer.
% 
% *Opts*: A list containing the optional input parameters. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.verbose: Flag for print out the number of bootstrap samples, 
% logical 0 or 1. Default value: 0.
% * Opts.rep: Number of replicates for scales. This option imposes special 
% structure on scaling parameters. For example, if Opts.rep = [3 4], this 
% means that the first three responses have the same scale and the next 
% four responses share a different scale. The elements of this vector should 
% sum to r. If not specified, the default is [], then all responses will be
% scaled differently. If all responses have the same scale, input [r], then 
% the regular envelope will be applied to the data.
% The input should be a row vector.
% * Opts.Gamma: The initial value for the envelope subspace. A p by u matrix. Default
% value is the one obtained from xenv, with $\Lambda^{-1}X$ being the 
% predictor and Y being the response. 
% * Opts.Lambda: The initial value for the scales. A p by p diagonal 
% matrix. Default value is the identity matrix. 
%
%% Output
%
% *bootse*: The standard error for elements in $\beta$ computed by
% bootstrap.  A p by r matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the scaled SIMPLS algorithm by bootstrapping the residuals. 

%% Example
%
%         load('chemo.mat')
%         X = X(:, [6 11 21 22]);
%         ModelOutput = spls(X, Y, 3);
% 
%         Opts.Gamma = ModelOutput.Gamma;
%         Opts.verbose = 1;
%         B = 10;
%         bootse = bstrp_spls(X, Y, 3, B, Opts)


function bootse = bstrp_spls(X, Y, u, B, Opts)

if nargin < 4
    error('Inputs: X, Y, u and B should be specified!');
elseif nargin == 4
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

X = double(X);
Y = double(Y);

[n, r] = size(Y);
p = size(X, 2);

ModelOutput = spls(X, Y, u, Opts);

Yfit = ones(n, 1) * ModelOutput.muY' + (X - ones(n, 1) * ModelOutput.muX') * ModelOutput.beta;
resi = Y - Yfit;
Opts.Lambda = ModelOutput.Lambda;
Opts.Gamma = ModelOutput.Gamma;
bootBeta = zeros(B, p * r);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of bootstrap sample ' int2str(i) '\n']);
    end
    
    bootresi = resi(randsample(1 : n, n, true), :);
    Yboot = Yfit + bootresi;
    temp = spls(X, Yboot, u, Opts);
    bootBeta(i, :) = reshape(temp.beta, 1, p * r);
    
end

bootse = reshape(sqrt(diag(cov(bootBeta, 1))), p, r);

fprintf('\nIf convergence is not reached for a bootstrap sample, \nit is still used in computing bootse.\n')
