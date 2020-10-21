function [tY, tL, rownan] = nk_RemNanSamples(Y, L, rownan)

if ~exist('rownan','var')
    
    nanmat = isnan(Y);
    rownan = any(nanmat,2);
    if sum(rownan)==0
        tY = Y;
        if exist('L','var'), tL = L; end
        return
    elseif sum(rownan) == size(Y,1)
        error('No samples in matrix after removal of NaN observations! Check your data or cross-validation structure')
    end
    tY = Y;
    tY(rownan,:) = [];
    if exist('L','var')
        tL = L;
        tL(rownan)=[];
    end
    
else
    
    tY = nan(numel(rownan),size(Y,2));
    tY(~rownan,:) = Y;
    if exist('L','var')
        tL = nan(numel(rownan),size(L,2));
        tL(~rownan,:) = L;
    end
    
end
