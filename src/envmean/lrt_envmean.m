%% lrt_envmean
% Select the dimension of the envelope subspace using likelihood ratio
% testing.

%% Syntax
%         u = lrt_envmean(Y, alpha)
%         u = lrt_envmean(Y, alpha, Opts)
%
%% Input
%
% *Y*: Data matrix. An n by p matrix, p is the dimension of the variable
% and n is number of observations. 
% 
% *alpha*: Significance level for testing.  A real number between 0 and 1,
% often taken at 0.05 or 0.01.
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
% This function implements the likelihood ratio testing procedure to select
% the dimension of the envelope subspace, with pre-specified significance 
% level $$\alpha$.  

%% Example
% 
%         load wheatprotein.txt
%         Y = wheatprotein(:, 1 : 6);
%         alpha = 0.01;
%         u = lrt_envmean(Y, alpha)

function u = lrt_envmean(Y, alpha, Opts)

if nargin < 2
    error('Inputs: Y and alpha should be specified!');
elseif nargin == 2
    Opts = [];
end

if (alpha < 0 || alpha > 1)
    error('alpha should be between [0, 1]!');
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n p] = size(Y);

ModelOutput0 = envmean(Y, p, Opts);

for i = 0 : p - 1

    	if printFlag == 1 
		fprintf(['Current dimension ' int2str(i) '\n']);
        end
    
	ModelOutput = envmean(Y, i, Opts);
	chisq = - 2 * (ModelOutput.l - ModelOutput0.l);
	df = ModelOutput0.np - ModelOutput.np;
	
    if chi2cdf(chisq, df) < (1 - alpha)
	    u = i;
	    break;
    end
    
end

if i == p - 1 && chi2cdf(chisq, df) > (1 - alpha)
    u = p;
end