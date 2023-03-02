close all;
clear all;
clc;
addpath('.\measure');
addpath('.\Database');
load('EYB_Group1.mat');
X=mapminmax(X,0,1);%0-1
gnd=labels;
K=max(gnd);
[d n]=size(X);


%% Parameters
lambda=1; 
p=0.95; 
c=0.11;
k=3;
alpha=4;

%% segmentation 

[Z,L,E,iter,EE] = solve_irfllrr(X,lambda,c,p,k);

%% postprocessing
[U,S,V] = svd(Z,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);
S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^(2*alpha);

%% spectral clustering NCut
idx = spectral_clustering(L, K);
Y=gnd;
predY=idx;
[result,bestY] = Clustering8Measure(Y, predY);
disp(['ACC nmi Purity Fscore Precision Recall AR Entropy=' num2str(result)]);

