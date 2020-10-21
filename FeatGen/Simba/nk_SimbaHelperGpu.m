function [pplus , pminus, delta_plus, delta_minus, marginList, dists] = ...
    nk_SimbaHelperGpu(g_X_train, g_X_train_square, g_Y_train,  Y_train, g_alphas, xs, labelsSet)

g_xs = GPUsingle(xs);
F = (g_alphas.^2)';

S = g_X_train(g_xs,:);
SL = g_Y_train(g_xs);
Pnum = size(g_X_train,1);
Snum = size(S,1);

% P2 =P.^2;
Pnorm = g_X_train_square * F; % sum(P.^2,2);
PNorm = Pnorm(:,1:Snum);
Snorm = ((S.^2) * F)'; % sum(S.^2,2)';
SNorm = zeros(Pnum,1,GPUsingle);
SNorm = SNorm+Snorm;

% Alphas = diag(alphas);
if Snum > 1,
    lenF = length(F);
    F = spdiags(F,0,lenF,lenF);
    SP = 2 * g_X_train * ( F * (S') );
else
    SP = 2 * g_X_train * (S'.* F);
end
g_dists = SNorm + PNorm - SP;
dists = single(g_dists);

for k = 1:Snum, 
    dists(xs(k),k) = realmax;  % to avoid P_pos to be the same point
end;

pplus = zeros(size(SL));
pminus = zeros(size(SL));
delta_plus = zeros(size(SL));
delta_minus = zeros(size(SL));

for Li = 1:length(labelsSet), 
    L                       = labelsSet(Li);
    Lcolumns                = (Y_train == L);
    FLcol                   = find(Lcolumns);
    FNLcol                  = find(~Lcolumns);
    Llines                  = find(SL == L);
    [min_pos, posI]         = min(dists(FLcol,Llines),[],1);
    pplus(Llines)           = FLcol(posI);
    min_pos                 = max(min_pos,0);
    delta_plus(Llines)      = sqrt(min_pos);
    [min_neg, negI]         = min(dists(FNLcol,Llines),[],1);
    pminus(Llines)          = FNLcol(negI);
    min_neg                 = max(min_neg,0);
    delta_minus(Llines)     = sqrt(min_neg);
end

marginList = delta_minus - delta_plus;
