function sY = apply_binning(Y, Bins)

[m,n] = size(Y);

sY = zeros(m,n);

for j=1:n
    vec = Bins{j};
    if isempty(vec), continue, end
    for k=1:numel(vec), sY( Y(:,j) >= vec(k) ) = k; end
end

end