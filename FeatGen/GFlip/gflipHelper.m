function [P_pos, P_neg, mindist_pos, mindist_neg, marginList, dists] = gflipHelper(X, Y, dists, changed_feat, added, cpus, labelsSet)

% [P_pos, P_neg, mindist_pos, mindist_neg, marginList, dists] = gflipHelper(X, X_square, Y,  alphas, inds, dists, add_feat, rem_feat);
%
% This function finds the nearhit and nearmiss for the instances in X(inds,:) with respect to X
% (while avoiding self nearhit).
%
% This function is used by simba in order o find the nearhit and nearmiss
% and the corresponding margins.
% 
% input: X(i,j) is the value of feature j in instance i
%        X_square = X.^2 (is given in order to save running time)
%        alphas(j) = the weights of the j's feature
%        inds = the indecies of the instances in interest
%
% output: P_pos = the indecies of the nearhits in X (i.e. X(P_pos(1)) is the nearhit for X(inds(1),:) )
%         P_neg = the indecies of the nearmiss in X (i.e. X(P_neg(1)) is the nearmiss for X(inds(1),:) )
%         mindist_pos(i) is the distance between X(inds(i),:) and its nearhit
%         mindist_neg(i) is the distance between X(inds(i),:) and its nearmiss
%         marginList(i) is the margin of X(inds(i),:), i.e. mindist_neg(i)-mindist_pos(i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Ran Gilad-Bachrach                    %%
%% Date: May 5, 2004                                              %%
%% Last update: May 6, 2004 by Amir Navot                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PL = Y;
%SL = PL;
m = size(X,1);
[y1, y2] = size(Y);
feat = X(:,changed_feat);
F1 = feat(:,ones(1,m));
F2 = F1';
delta = (F1 - F2).^2;
if added,
    dists = dists + delta;
else
    dists = dists - delta;
end;

% dists = max(dists,diag(realmax*ones(1,m))); % to avoid choosing an instance as the closest prototype to itself

%labelsSet = unique(PL);

P_pos = zeros(y1,y2);
P_neg = zeros(y1,y2);
mindist_pos = zeros(y1,y2);
mindist_neg = zeros(y1,y2);

for Li = 1:length(labelsSet), 
  
    L = labelsSet(Li);
    Lcolumns = (Y == L);
    FLcol = find(Lcolumns);
    FNLcol = find(~Lcolumns);
    %Llines = find(Y == L);
  
    [min_pos, posI] = gflipmin(dists(FLcol,Lcolumns),1,cpus);
    %[min_pos, posI] = min(dists(FLcol,Llines),[],1);
    P_pos(Lcolumns) = FLcol(posI);
    min_pos = max(min_pos,0);
    mindist_pos(Lcolumns) = sqrt(min_pos);
    
    [min_neg, negI] = gflipmin(dists(FNLcol,Lcolumns),1,cpus);
    %[min_neg, negI] = min(dists(FNLcol,Llines),[],1);
    P_neg(Lcolumns) = FNLcol(negI);
    min_neg = max(min_neg,0);
    mindist_neg(Lcolumns) = sqrt(min_neg);
end

marginList = mindist_neg - mindist_pos;
%suggestedBeta = 1 / mean(marginList);
%[margin] = min(marginList);






