function bootse = bstrp_envmean(X, u, B, Opts)

if nargin < 3
    error('Inputs: X, B and u should be specified!');
elseif nargin == 3
    Opts = [];
end

Opts = make_opts(Opts);
printFlag = Opts.verbose;
Opts.verbose = 0;

[n p] = size(X);

ModelOutput = envmean(X, u, Opts);

Xmean = ones(n, 1) * ModelOutput.mu';
resi = X - Xmean;

bootBeta = zeros(B, p);

for i = 1 : B
    
    if printFlag == 1
        fprintf(['Current number of bootstrap sample ' int2str(i) '\n']);
    end
    
    bootresi = resi(randsample(1 : n, n, true), :);
    Xboot = Xmean + bootresi;
    temp = envmean(Xboot, u, Opts);
    bootBeta(i, :) = temp.mu';

end

bootse = sqrt(diag(cov(bootBeta, 1)));