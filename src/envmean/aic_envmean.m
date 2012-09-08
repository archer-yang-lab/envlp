%% aic_envmean
% Select the dimension of the envelope subspace using Akaike information
% criterion.

%% Syntax
%         u = aic_envmean(Y)
%         u = aic_envmean(Y, Opts)
%
%% Input
%
% *Y*: Data matrix. An n by p matrix, p is the dimension of the variable
% and n is number of observations. 
% 
% *Opts*: A list containing the optional input parameters, to control the
% iterations in sg_min. If one or several (even all) fields are not
% defined, the default settings are used.
% 
% * Opts.maxIter: Maximum number of iterations.  Default value: 300.
% * Opts.ftol: Tolerance parameter for F.  Default value: 1e-10. 
% * Opts.gradtol: Tolerance parameter for dF.  Default value: 1e-7.
% * Opts.verbose: Flag for print out dimension selection process, 
% logical 0 or 1. Default value: 0.
%
%% Output
%
% *u*: Dimension of the envelope. An integer between 0 and p.

%% Description
% This function implements the Akaike information criteria (AIC) to select
% the dimension of the envelope subspace. 

%% Example
%         load wheatprotein.txt
%         Y = wheatprotein(:, 1 : 6);
%         u = aic_envmean(Y)

function u = aic_envmean(Y, Opts)

if nargin < 1
    error('Inputs: Y should be specified!');
elseif nargin == 1
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n p] = size(Y);
    
ModelOutput = envmean(Y, p, Opts);
ic = - 2 * ModelOutput.l + 2 * ModelOutput.paramNum;
u = p;

for i = 0 : p - 1
    
	if printFlag == 1 
		fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
	ModelOutput = envmean(Y, i, Opts);
	temp = -2 * ModelOutput.l + 2 * ModelOutput.paramNum;
	
    if temp < ic
	   u = i;
	   ic = temp;
    end
    
end
