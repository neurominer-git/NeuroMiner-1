function W = nk_WeigthDataInstanceHisto(L)

% %% Compute target label histogram
% nbin = hist(L);
% cbin = min(L) : ( max(L) - min(L) ) / numel(nbin) : max(L);
% N = numel(L);
% 
% %% Compute inverse relative histogram for weighting
% nbin = nk_ScaleData(1 - nbin / N,1,10); 
% W = zeros(N,1);
% 
% %% Assign weights according to bin info
% for i = 1:numel(cbin)-1
%    ind = L >= cbin(i) & L < cbin(i+1);
%    W(ind) = nbin(i);
% end
% ind = L == cbin(i+1);
% W(ind) = nbin(i);

W = nk_ScaleData(abs(L - median(L)),0,1);