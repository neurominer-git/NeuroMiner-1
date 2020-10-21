function G = nk_CompGridFreq(V1, V1ind, V2, V2ind)

vl1 = length(V1ind);
vl2 = length(V2ind);
G = zeros(vl1,vl2);
G1 = zeros(length(V1),1);
G2 = zeros(length(V2),1);

for z=1:length(V1ind)
    ind = V1 == V1ind(z);
    G1(ind) = z;
end
for z=1:length(V2ind)
    ind = V2 == V2ind(z);
    G2(ind) = z;
end

for z=1:length(V1)
    try
        G(G1(z),G2(z)) =  G(G1(z),G2(z)) +1;
    catch
         G(1,G2(z)) =  G(1,G2(z)) +1;
    end
end