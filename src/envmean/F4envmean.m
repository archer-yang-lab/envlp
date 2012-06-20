function f = F4envmean(R, DataParameter)

n = DataParameter.n;
p = DataParameter.p;
sX = DataParameter.sX;
sigX = DataParameter.sigX;
logDetSX = DataParameter.logDetSX;

eigtem = eig(R' * sigX * R);
a = log(prod(eigtem(eigtem > 0)));

eigtem0 = eig(R' * inv(sX) * R);
b = log(prod(eigtem0(eigtem0 > 0)));

f = n * p * (1 + log(2 * pi)) + n * (a + b + logDetSX);