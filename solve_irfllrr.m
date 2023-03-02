function [Z,L,E,iter,EE] = solve_irfllrr(X,lambda,c,p,k)
if nargin<2
    lambda = 1.0; 
end
%% 
Q1 = orth(X');
Q2 = orth(X);
A = X*Q1;
B = Q2'*X;

[Z,L,E,iter,EE] = solve_irfllrra(X,A,B,lambda,c,p,k);
Z = Q1*Z;
L = L*Q2';