function [ Lperm, Yperm ] = nk_VisXPermY(Y, L, IND, permmode, indpermrows, indpermcols, permind)
global MODEFL

switch permmode
    
    case 1
        Lperm = L(indpermrows(:,permind));
        indN = isnan(Lperm); 
        uL = unique(Lperm(~indN)); nuL = numel(uL);
        Yperm = Y;
    case {2,3}
        if permmode == 3, 
            Lperm = L(indpermrows(:,permind)); 
        else
            Lperm = L(IND);
        end
        L = L(IND);
        indN = isnan(L);
        uL = unique(L(~indN)); nuL = numel(uL);
        Yperm = zeros(size(Y));
        for i=1:nuL
           indi = L == uL(i);
           Yperm(indi,:) = Y(indi,indpermcols(i,:,permind));
        end
end

if strcmp(MODEFL,'classification')
    tL = Lperm;
    for i=1:nuL
        if i==1,
            tL( Lperm == uL(1) ) = 1;
        else
            tL( Lperm == uL(i) ) = -1;
        end
    end
    Lperm = tL;
end
if sum(indN)
   Lperm(indN)=[]; Yperm(indN,:)=[]; 
end
