function k=Kpd(p,d)

k=zeros(p*d,p*d);

for i=1:p
    for j=1:d
        Hij=zeros(p,d);
        Hij(i,j)=1;
        k=k+kron(Hij,Hij');
    end
end