clear
load('/home/yifeng/YifengLi/Research/ML/srv1_7/data/ColonCancer.mat');
dataStr='Colon';
classes(classes==0)=-1;
X=D;
L=classes;
k=8;
[AtA,Y,numIter,tElapsed,finalResidual]=superkernelseminmfnnls(X,L,k);

