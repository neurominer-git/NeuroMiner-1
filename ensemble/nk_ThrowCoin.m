function H = nk_ThrowCoin(H)
hi = H==0;
if any(hi) > 0
    hx0 = hi; 
    shx0 = sum(hx0);
    i0 = rand(shx0,1);
    i0(i0>0.5) = 1; i0(i0~=1)=-1;
    H(hx0) = i0;
end
return