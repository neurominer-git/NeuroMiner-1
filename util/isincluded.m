function v = isincluded(Subset, Set)
% ISINCLUDED (Subset, Set) returns a vector v of length(Set) 
% logical values. A 1 is returned for every element of Set that is in
% Subset, 0 otherwise.
% 
% Example isincluded([1 2], [0 1 1 2 4 5 1 2 4 5]) returns
% [0 1 1 1 0 0 1 1 0 0]

v = zeros (length(Set),1);
for i=1:length(Subset)
    ii = find(Subset(i) == Set);
    v(ii) = 1;
end
v=logical(v); 