function [C, Param] = nk_SortComponents(Param, vecname, actind, shelfind)

[ ix, jx ] = size(Param);
n = ix*jx; ll = 2;
mp = Param{1}{actind}{shelfind}.mpp.(vecname);
while ll<=n
    mpll = Param{ll}{actind}{shelfind}.mpp.(vecname);
    C = zeros(size(mp,2),size(mpll,2));
    for j = 1:size(C,2)
        C(:,j) = nk_CorrMat(mp, mpll(:,j));
    end
    ll=ll+1;
end

