function u=bic_env(Y,X)

[n r]=size(Y);
p=size(X,2);
    
[beta Sigma Gamma Gamma0 eta Omega Omega0 alpha l ratio]=env(X,Y,r);
ic=-2*l+log(n)*(r*p+r*(r+1)/2);
u=r;


for i=0:r-1
    i
        [beta Sigma Gamma Gamma0 eta Omega Omega0 alpha l ratio]=env(X,Y,i);
        temp=-2*l+log(n)*(u*p+r*(r+1)/2);
        if (temp<ic)
           u=i;
           ic=temp;
        end
end
