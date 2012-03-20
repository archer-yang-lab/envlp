function k=Kpd_right(X,p,d)

% p is the number of row and d is number of column before vectorize the
% matrix
n=size(X,1);
m=size(X,2);
k=zeros(size(X));
for i=1:n
    temp=reshape(X(i,:),d,p);
    temp=temp';
    k(i,:)=reshape(temp,1,m);
end