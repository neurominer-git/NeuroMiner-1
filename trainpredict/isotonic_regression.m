function ghat = isotonic_regression(y,w)
% ghat = IsoMeans(y)
% fits a vector ghat with nondecreasing components to the 
% data vector y such that sum((y - ghat).^2) is minimal. 
% (Pool-adjacent-violators algorithm).
% In case of a second input vector w with positive entries 
% (and the same size as y),
% ghat = IsoMeans(y,w)
% produces an isotonic vector minimizing sum(w.*(y - ghat).^2) 
% 
% Lutz Duembgen, 23.02.2000

n = length(y);
if nargin == 1
	w = ones(size(y));
end
index  = zeros(size(y));
weight = zeros(size(y));
% An interval of indices is represented by its left endpoint 
% ("index") and its "weight" (sum of w's).
ghat = zeros(size(y));

ci = 1;
index(ci) = 1;
weight(ci) = w(1);
ghat(ci) = y(1);
% ci is the number of the interval considered currently.
% ghat(ci) is the weighted mean of y-values within this interval.
for j=2:n
	% a new index intervall, {j}, is created:
	ci = ci+1;
	index(ci) = j;
	weight(ci) = w(j);
	ghat(ci) = y(j);
	while ci >= 2 && ghat(ci-1) >= ghat(ci)
		% "pool adjacent violators":
		nw = weight(ci-1) + weight(ci);
		ghat(ci-1) = (weight(ci-1)*ghat(ci-1) + weight(ci)*ghat(ci)) ...
			/ nw;
		weight(ci-1) = nw;
		ci = ci-1;
	end
end
% Now define ghat for all indices:
while n >= 1
	for j=index(ci):n
		ghat(j) = ghat(ci);
	end;
	n = index(ci)-1;
	ci = ci-1;
end
end