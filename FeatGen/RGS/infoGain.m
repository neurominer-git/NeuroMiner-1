function scores = infoGain(X_train, Y_train)

% scores = infoGain(X_train, Y_train);
%
% Implementation of infoGain feature selection, as used in:
%
% A. Navot, L. Shpigelman, N. Tishby, E. Vaadia. Nearest Neighbor Based Feature Selection for Regression and its
% Application to Neural Activity. Submitted to NIPS 2005.
%
% input: X_train(i,j) is the value of feature j in training instance i.
%        Y_train(i) is the label of training instance i.
%
% output: scores(j) is the mutual information between the j'th feature and the target values. higher is better.
%                  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Lavi Shpigelman               
%% Date: June 3, 2005                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% force the data to be binary
Y_train = Y_train > median(Y_train);
med = median(X_train);
X_train = X_train > med(ones(1,size(X_train,1)),:);

scores = featureInformation(X_train, Y_train);


function fMI = featureInformation(X,y);
feat_num = size(X,2);
X = X > 0;
y = y > 0;
y = y(:);
fMI = zeros(1,feat_num);
for fi = 1:feat_num,
    t = 2*X(:,fi)+y;
    p = [sum(t==0), sum(t==1); sum(t==2), sum(t==3)];
    p = p / length(t);
    fMI(fi) = MutualInformation(p);
end

function MI = MutualInformation(pxy);
px = sum(pxy,2);
py = sum(pxy);
pxpy = px * py;
MI = dkl(pxy(:), pxpy(:));

function d = dkl(p, q);
logpq = zeros(size(p));
I = find(p>eps);
logpq(I) = log2(p(I)./q(I));
d = (sum(p .* logpq));