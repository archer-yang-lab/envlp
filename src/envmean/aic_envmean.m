function u = aic_envmean(X, Opts)

if nargin < 1
    error('Inputs: X should be specified!');
elseif nargin == 1
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n p] = size(X);
    
ModelOutput = envmean(X, p, Opts);
ic = - 2 * ModelOutput.l + 2 * ModelOutput.np;
u = p;

for i = 0 : p - 1
    
	if printFlag == 1 
		fprintf(['Current dimension ' int2str(i) '\n']);
    end
    
	ModelOutput = envmean(X, i, Opts);
	temp = -2 * ModelOutput.l + 2 * ModelOutput.np;
	
    if temp < ic
	   u = i;
	   ic = temp;
    end
    
end
