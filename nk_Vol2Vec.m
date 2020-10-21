function Y = nk_Vol2Vec(P)

if ~exist('P','var') || ~exist(P,'file')
    
    P = spm_select(Inf, 'image', 'Select images');
    
end

V = spm_vol(P);
Y = zeros(numel(V),prod(V(1).dim));

for i = numel(V);
   y = spm_read_vols(V(i)); 
   Y(i,:) = reshape(y,1,numel(y));
end

end