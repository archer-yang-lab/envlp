%% bstrp_henv
% Compute bootstrap standard error for the heteroscedastic envelope model. 

%% Usage
% bootse=bstrp_henv(X,Y,B,u)
%
% Input
%
% * X: Group indicators. An n by p matrix, p is the number of groups. X can
% only take p different values, one for each group.
% * Y: Multivariate responses, an n by r matrix, r is the number of
% responses and n is number of observations.  The responses must be continuous variables.
% * B: Number of boostrap samples.  A positive integer.
% * u: Dimension of the envelope subspace.  A positive integer between 0 and
% r.
%
% Output
%
% * bootse: The standard error for elements in $$\beta$ computed by
% bootstrap.  An r by p matrix.

%% Description
% This function computes the bootstrap standard errors for the regression
% coefficients in the heteroscedastic envelope model by bootstrapping the residuals. 

function bootse=bstrp_henv(X,Y,B,u)

dataParameter=make_parameter4henv(X,Y);
p=dataParameter.p;
r=dataParameter.r;
n=dataParameter.n;
ng=dataParameter.ng;
ncum=dataParameter.ncum;
ind=dataParameter.ind;

stat=henv(X,Y,u);

Yfit=stat.Yfit;
resi=Y-Yfit;

bootresi=zeros(n,r);
bootBeta=zeros(B,r*p);

for i=1:B
    
    for j=1:p
        if j>1
            bootresi(ind(ncum(j-1)+1:ncum(j)),:)=resi(randsample(ind(ncum(j-1)+1:ncum(j)),ng(j),true),:);   
        else
            bootresi(ind(1:ncum(1)),:)=resi(randsample(ind(1:ng(1)),ng(1),true),:);
        end
    end
    Yboot=Yfit+bootresi;
    temp=henv(X,Yboot,u);
    bootBeta(i,:)=reshape(temp.beta,1,r*p);
    
end

bootse=reshape(sqrt(diag(cov(bootBeta,1))),r,p);