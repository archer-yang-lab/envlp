function u=lrt_env(Y,X,alpha)

[n r]=size(Y);
p=size(X,2);

[beta Sigma Gamma Gamma0 eta Omega Omega0 alpha f0 ratio]=env(X,Y,r);


for i=0:r
    i
        [beta Sigma Gamma Gamma0 eta Omega Omega0 alpha f ratio]=env(X,Y,i);
        chic = -2*(f-f0);
        
        if (chi2cdf(chic,p*(r-i)) < (1-alpha))
            u=i;
            break;
        end
end

