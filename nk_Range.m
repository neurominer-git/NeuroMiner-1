function R = nk_Range(X, dim)

if ~exist('dim','var'), dim=1; end

min_X = min(X,[],dim);
max_X = max(X,[],dim);
R = abs(max_X - min_X);