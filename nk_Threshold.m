function [T, Thresh] = nk_Threshold(V, thresh)

if thresh.type ~= 3
    sig = false(size(V,2),1); 
    switch thresh.type
        case 1
            % Percentile method
            Thresh = percentile(V, thresh.val);
        case 2
            % Absolute threshold
            Thresh = thresh.val;
           
    end
    switch thresh.logop
        case 1
            ind = V > Thresh;
        case 2
            ind = V >= Thresh;
        case 3
            ind = V < Thresh;
        case 4
            ind = V <= Thresh;
        case 5
            ind(:,1) = V < max(Thresh);
            ind(:,2) = V > min(Thresh);
        case 6
            ind(:,1) = V <= max(Thresh);
            ind(:,2) = V >= min(Thresh);
        case 7
            ind(:,1) = V > max(Thresh);
            ind(:,2) = V < min(Thresh);
        case 8
            ind(:,1) = V >= max(Thresh);
            ind(:,2) = V <= min(Thresh);

    end
    for i=1:size(ind,2), sig(ind(:,i)) = true; end
    T = zeros(size(V));
    T(sig) = V(sig);
else
    T=V;Thresh = [];
end


return