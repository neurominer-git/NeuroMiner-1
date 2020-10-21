% examples of using universe sparse representation (versatile sparse matrix factorization: VSMF)

clear
% USR/VSMF
load('/home/yifeng/YifengLi/Research/ML/srv1_7/data/smallRoundBlueCellTumorOFChildhood.mat');
dataStr='SRBCT';
D=Data;
clear('Data');

classes=changeClassLabels01(classes);

optionusr.alpha2=0;
optionusr.alpha1=0.01;
optionusr.lambda2=0;
optionusr.lambda1=0.01;
optionusr.t1=false;
optionusr.t2=true;
k=6;
% [A,Y]=sparsenmf2rule(D,k,optionsnmf2);
[A,Y]=vsmf(D,k,optionusr);
