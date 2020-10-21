function suggestedBeta = suggestBeta(X, Y);

% suggestedBeta = suggestBeta(obj1, X, Y);
%
% This function find a reasonable value of beta for the sigmoid utility.
% You can use this value as starting point in tuning this parameter.
%
% input: X(i,j) is the value of feature j in instance i
%        Y(i) is the label of instance i
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot                                         %%
%% Date: May 8, 2005                                              %%
%% Last update: May 8, 2005 by Amir Navot                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = X;
P2 = X.^2;
PL = Y;

inds = 1:length(Y);
alphas = ones(size(X,2),1);

S = P(inds,:);
SL = PL(inds);

Pnum = size(P,1);
Snum = size(S,1);

% P2 =P.^2;
Pnorm = P2 * alphas; % sum(P.^2,2);
PNorm = Pnorm(:,ones(1, Snum));
Snorm = ((S.^2) * alphas)'; % sum(S.^2,2)';
SNorm = Snorm(ones(1,Pnum), :);

% Alphas = diag(alphas);
if Snum > 1,
    Alphas = spdiags(alphas,0,length(alphas),length(alphas));
    SP = 2 * P * ( Alphas * (S') );
else,
    SP = 2 * P * (S'.*alphas);
end

dists = SNorm + PNorm - SP;

for k = 1:Snum, 
    dists(inds(k),k) = realmax;  % to avoid P_pos to be the same point
end;

if any(dists(:) == 0),
    warning('There are identical instances in the data set');
end

labelsSet = unique(PL);

P_pos = zeros(size(SL));
P_neg = zeros(size(SL));
mindist_pos = zeros(size(SL));
mindist_neg = zeros(size(SL));


for Li = 1:length(labelsSet), 
  L = labelsSet(Li);
  Lcolumns = (PL == L);
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
suggestedBeta = 1 / mean(marginList);







