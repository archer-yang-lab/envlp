




function u=lrt_ienv(X,Y,alpha)

[n r]=size(Y);
p=size(X,2);

stat0=env(X,Y,r);


for i=1:(p+1)
%     i
        stat=ienv(X,Y,p+1-i);
        chisq = -2*(stat.l-stat0.l);
        df=stat0.np-stat.np;
        
        if (chi2cdf(chisq,df) < (1-alpha))
            u=p+1-i;
            break;
        end
end

if (i== p+1) && chi2cdf(chisq,df) > (1-alpha)
    u=r;
    warning('No inner envelope model is selected, fit with the standard multivariate linear model.');
end