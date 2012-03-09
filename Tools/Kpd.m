function k=Kpd(X,p,d)

% p is the number of row and d is number of column before vectorize the
% matrix
n=size(X,1);
m=size(X,2);
k=zeros(size(X));
for i=1:m
    temp=reshape(X(:,i),p,d);
    temp=temp';
    k(:,i)=reshape(temp,n,1);
end