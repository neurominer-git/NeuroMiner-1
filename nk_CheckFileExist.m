function [PP, fPP] = nk_CheckFileExist(P)

PP=[]; fPP=[];

if iscell(P), P = char(P);end

for i=1:size(P,1)
    iP = deblank(P(i,:));
    if ~exist(iP,'file')
        if isempty(PP)
            PP = iP;
        else
            PP = char(PP,iP);
        end
    else
        if isempty(fPP)
            fPP = iP;
        else
            fPP = char(fPP,iP);
        end
    end
    
end