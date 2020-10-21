function [P_pos, P_neg, mindist_pos, mindist_neg, marginList, dists] = nk_simbaHelper(X, X_square, Y,  alphas, inds, labelsSet)

% [P_pos, P_neg, mindist_pos, mindist_neg, marginList, dists] = simbaHelper(X, X_square, Y,  alphas, inds);
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
%         dists(i1,i2) is the square distance between X(inds(i1),:) and X(i2,:)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Ran Gilad-Bachrach                    %%
%% Date: April 1, 2004                                            %%
%% Last update: May 6, 2004 by Amir Navot &    Ran Gilad-Bachrach %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = X(inds,:);
SL = Y(inds);

alphas = alphas(:);
Pnum = size(X,1);
Snum = size(S,1);

% P2 =P.^2;
Pnorm = X_square * alphas; % sum(P.^2,2);
PNorm = Pnorm(:,ones(1, Snum));
Snorm = ((S.^2) * alphas)'; % sum(S.^2,2)';
SNorm = Snorm(ones(1,Pnum), :);

% Alphas = diag(alphas);
if Snum > 1,
    Alphas = spdiags(alphas,0,length(alphas),length(alphas));
    SP = 2 * X * ( Alphas * (S') );
else
    SP = 2 * X * (S'.* alphas);
end

dists = SNorm + PNorm - SP;

for k = 1:Snum, 
    dists(inds(k),k) = realmax;  % to avoid P_pos to be the same point
end;

P_pos = zeros(size(SL));
P_neg = zeros(size(SL));
mindist_pos = zeros(size(SL));
mindist_neg = zeros(size(SL));
for Li = 1:length(labelsSet), 
  L = labelsSet(Li);
  Lcolumns = (Y == L);
  FLcol = find(Lcolumns);
  FNLcol = find(~Lcolumns);
  Llines = find(SL == L);
  [min_pos, posI] = min(dists(FLcol,Llines),[],1);
  P_pos(Llines) = FLcol(posI);
  min_pos=max(min_pos,0);
  mindist_pos(Llines) = sqrt(min_pos);
  [min_neg, negI] = min(dists(FNLcol,Llines),[],1);
  P_neg(Llines) = FNLcol(negI);
  min_neg=max(min_neg,0);
  mindist_neg(Llines) = sqrt(min_neg);
end

marginList = mindist_neg - mindist_pos;
%[margin] = min(marginList);






