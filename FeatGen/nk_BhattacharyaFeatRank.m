function D = nk_BhattacharyaFeatRank(Y, L)

% Y     :       Data
% L     :       Target Labels

[~,n] = size(Y);
ix = unique(L(~isnan(L)));
indP = L==ix(1); indM = L==ix(2);
YP = Y(indP,:); YM = Y(indM,:);
D = zeros(1,n);
warning off
for i=1:size(Y,2)
   ind = true(1,n); ind(i)=false;
   D(i) = bhattacharyya(YP(:,ind),YM(:,ind));
end
D=1./D';
