function u=bic_ienv(X,Y)

[n r]=size(Y);
p=size(X,2);
    
stat=env(X,Y,r);
ic=-2*stat.l+log(n)*stat.np;
u=r;


for i=0:p
%     i
        stat=ienv(X,Y,i);
        temp=-2*stat.l+log(n)*stat.np;
        if (temp<ic)
           u=i;
           ic=temp;
        end
end
