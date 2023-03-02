function [Z,L,E,iter,EE] = solve_irfllrra(X,A,B,lambda,c,p,k)
rho=1.12;
tol = 1e-6;   
eps=1e-16;
maxIter = 1e6;
[d n] = size(X);
nA = size(A,2);
dB = size(B,1);
max_mu = 1e10;
mu = 1e-6;
ata = A'*A;
bbt = B*B';

%% Initializing optimization variables
Z = zeros(nA,n);
L = zeros(d,dB);
E = sparse(d,n);

Y1 = zeros(d,n);
Y2 = zeros(nA,n);
Y3 = zeros(d,dB);
%% Start main loop
iter = 0;
display=1;
disp(['initial,r(Z)=' num2str(rank(Z)) ',r(L)=' num2str(rank(L)) ',|E|_1=' num2str(sum(sum(abs(E))))]);

while iter<maxIter
    iter = iter + 1;
   
    %Update J
    temp1=Z+Y2/mu;
    [U,sigma,V] = svd(full(temp1),'econ');
    dsigma=diag(sigma);
    [sigma,svp]=reweighted(X,dsigma,c,p,mu,k,eps);
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    
    %Update S
    temp2=L+Y3/mu;
    [U,sigma,V] = svd(full(temp2),'econ');
    dsigma=diag(sigma);
    [sigma,svp]=reweighted(X,dsigma,c,p,mu,k,eps);
    S = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
     
    %udpate Z
    Z = inv(ata+eye(nA))*(A'*(X-L*B-E)+J+(A'*Y1-Y2)/mu);
    
    %update L
    L = ((X-A*Z-E)*B'+S+(Y1*B'-Y3)/mu)*inv(bbt+eye(dB));
    
    %update E
    temp3 = X-A*Z-L*B+Y1/mu;
    E = max(0,temp3 - lambda/mu)+min(0,temp3 + lambda/mu); %L-1
%     E = solve_l1l2(temp3,lambda/mu); %L-2

    %update the multiplies
    leq1 = X-A*Z-L*B-E;
    leq2 = Z-J;
    leq3 = L-S;
    stopC =max(max(max(abs(leq3))),max(max(max(abs(leq1))),max(max(abs(leq2)))));
    if display&&(iter==1 || mod(iter,50)==0 || stopC<tol)    
        disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ',stopALM=' num2str(stopC,'%2.3e') ...
            ',|E|_1=' num2str(sum(sum(abs(E))))]);
    end
    if stopC<tol
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        Y3 = Y3 + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
    EE(1,iter)=stopC;
end
function [E] = solve_l1l2(W,lambda)
n = size(W,2);%W的列数
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end

function [x] = solve_l2(w,lambda)%lemma4.1
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);%返回2范式
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end