function [ sigma,svp ] = reweighted(X,sigma,c,p,mu,k,eps)
[d,n]=size(X);
c=sqrt(d*n)*c;
temp=sigma;
for i=1:k
    W=c*(temp+eps).^(p-2);
    W=1./(W./mu+1);
    sigma=W.*sigma;
    sigma=max(sigma-eps,0);
    temp=sigma;
end
svp=length(find(sigma>0));
sigma=sigma(1:svp);
end
