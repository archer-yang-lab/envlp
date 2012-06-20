function df = dF4envmean(R, DataParameter)

n = DataParameter.n;
sX = DataParameter.sX;
sigX = DataParameter.sigX;

a = 2 * sigX * R * inv(R' * sigX * R);

temp = inv(sX);

b = 2 * temp * R * inv(R' * temp * R);

df = n * (a + b);