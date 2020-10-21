function scores = corrcoef_select(X_train, Y_train);  

% scores = corrcoef_select(X_train, Y_train);
%
% Implementation of correletion coeficient feature selection, as used in:
%
% A. Navot, L. Shpigelman, N. Tishby, E. Vaadia. Nearest Neighbor Based Feature Selection for Regression and its
% Application to Neural Activity. Submitted to NIPS 2005.
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%
% output: scores(j) is the squared correletaion coefficient between the j'th feature and the target values. higher is better.
%                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Lavi Shpigelman               
%% Date: June 3, 2005                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std1=std(Y_train,1);
std2=std(X_train,1);
num=(Y_train-mean(Y_train))'*(X_train-ones(size(X_train,1),1)*mean(X_train));
I=std2~=0;
W=zeros(1,size(X_train,2));
W(I)=num(I)./(std1*std2(I)*size(X_train,1));

scores=W.^2;