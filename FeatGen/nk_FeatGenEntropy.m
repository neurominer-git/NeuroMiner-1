function H = nk_FeatGenEntropy(X, dim, skipzeros)

if dim==2;
    X=X';
end

[Cl, kFea] = size(X); % # features , # observations
H = zeros(kFea,1);

for q=1:kFea
    if skipzeros 
        Xx = X(X(:,q)~=0,q);
        Cl = length(Xx);
    else
        Xx = X(:,q);
    end
    Cu = unique(Xx); 
    Cux = length(Cu); 
    Hx = 0; 
    for r=1:Cux
        P = sum(Xx == Cu(r))/Cl; 
        Hx = Hx + P*log2(P);
    end
    H(q) = -1*Hx;
end

