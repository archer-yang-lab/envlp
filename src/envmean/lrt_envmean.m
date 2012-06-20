function u = lrt_envmean(X, alpha, Opts)

if nargin < 2
    error('Inputs: X and alpha should be specified!');
elseif nargin == 2
    Opts = [];
end

if (alpha < 0 || alpha > 1)
    error('alpha should be between [0, 1]!');
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n p] = size(X);

ModelOutput0 = envmean(X, p, Opts);

for i = 0 : p - 1

    	if printFlag == 1 
		fprintf(['Current dimension ' int2str(i) '\n']);
        end
    
	ModelOutput = envmean(X, i, Opts);
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