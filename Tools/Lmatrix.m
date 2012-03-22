function L=Lmatrix(r)

L=zeros(r^2,r-1);
for i=1:r-1
    L(i*r+i+1,i)=1;
end