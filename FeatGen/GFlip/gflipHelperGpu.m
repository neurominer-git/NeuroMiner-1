function [P_pos, P_neg, mindist_pos, mindist_neg, marginList, g_dists] = gflipHelperGpu(X, Y, g_dists, changed_feat, added, labelsSet)

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


m = size(X,1);
feat = X(:,changed_feat);
F1 = feat(:,ones(1,m));
delta = F1-F1';
delta = delta.*delta;
if added,
    g_dists = g_dists + delta;
else
    g_dists = g_dists - delta;
end;

% dists = max(dists,diag(realmax*ones(1,m))); % to avoid choosing an
% instance as the closest prototype to itself

P_pos = GPUsingle(zeros(size(Y)));
P_neg = GPUsingle(zeros(size(Y)));
mindist_pos = GPUsingle(zeros(size(Y)));
mindist_neg = GPUsingle(zeros(size(Y)));
indvec = 1:size(Y,1);

for Li = 1:length(labelsSet), 
  
    L = labelsSet(Li);
    Lcolumns = (Y == L);
    FLcol = GPUsingle(indvec(Lcolumns));
    FNLcol = GPUsingle(indvec(~Lcolumns));
    g_distsFLcol = g_dists(FLcol,FLcol);
    g_distsFNLcol = g_dists(FNLcol,FLcol);
    [min_pos, posI] = gmin(g_distsFLcol,1);
    [min_pos, posI] = min(dists(FLcol,Lcolumns),[],1);
    P_pos(Lcolumns) = FLcol(posI);
    min_pos = max(min_pos,0);
    mindist_pos(Lcolumns) = sqrt(min_pos);
    
    [min_neg, negI] = min(dists(FNLcol,Lcolumns),[],1);
    P_neg(Lcolumns) = FNLcol(negI);
    min_neg = max(min_neg,0);
    mindist_neg(Lcolumns) = sqrt(min_neg);
end
marginList = mindist_neg - mindist_pos;
%suggestedBeta = 1 / mean(marginList);
%[margin] = min(marginList);

function [minvals, ind] = gmin(arr, dim)

switch dim
    case 1
        ind = GPUsingle(zeros(size(arr,2),1));
        N1 = size(arr,2);
        N2 = size(arr,1);
        for l=1:N1
            ind(l) = cublasIsamin(N2, getPtr(arr(:,l)), 1) + N2*(l-1);
        end
    case 2
        ind = GPUsingle(zeros(size(arr,1),1));
        N1 = size(arr,1);
        N2 = size(arr,2);
        for l=1:N1
            ind(l) = cublasIsamin(N2, getPtr(arr(l,:)), 1) + N2*(l-1);
        end
end

minvals = arr(ind);
   