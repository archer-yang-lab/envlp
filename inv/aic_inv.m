function u=aic_inv(X,Y)

[n r]=size(Y);
    
stat=env(X,Y,r);
ic=-2*stat.l+2*stat.np;
u=r;


for i=0:p
%     i
        stat=inv(X,Y,i);
        temp=-2*stat.l+2*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end