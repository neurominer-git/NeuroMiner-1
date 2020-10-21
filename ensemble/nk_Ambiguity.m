function A = nk_Ambiguity(P,~)
%
% Compute the entropy-based ambiguity according to the formula given in 
% Padraig Cunningham TCD-CS-2000-07
m = size(P,1);

Mx = 0;
for i=1:m % Loop through samples x_i
    Kx = 0;
    pos = P(i,:)~=0;    
    minx = min(P(i,pos));
    maxx = max(P(i,pos));
    for j=minx:maxx
        % Compute frequency of k_j (class) for sample x_i
        f = sum(P(i,pos)==j)/sum(pos);
        if ~f, continue, end
        Kx = Kx + (-f*log(f));
    end
    Mx = Mx + Kx;
end
A = Mx/m;
