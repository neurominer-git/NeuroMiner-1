function [Y, label, I, NanStats] = nk_ManageNanCases(Y, label, I, act)
% =========================================================================
% FORMAT [Y, label, I] = nk_ManageNanCases(Y, label, I)
% =========================================================================
% Removes cases that consist only of NaN from the data/label matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 09/2017
if ~exist('label','var'), label = []; end
if ~exist('act','var') || isempty(act) , act = 'prune'; end
switch act
    
    case 'prune'
        NanStats = [];
        if ~exist('I','var') || isempty(I) 
            % Remove NaN cases
            [~,n] = size(Y); 
            I = sum(isnan(Y),2) == n ; 
            if any(I), 
                Y(I,:) = []; 
                if ~isempty(label), 
                    label(I,:)=[]; 
                end; 
            end
        elseif ~isempty(I) && any(I)
            % Add NaN cases
            if ~isempty(Y),
                if iscell(Y),
                    for i=1:numel(Y), Y{i} = FillWithNan(Y{i}, I); end; 
                else
                    Y = FillWithNan(Y, I);
                end
            end
            if ~isempty(label); label = FillWithNan(label, I); end
        end
        
    case 'inform'
        [m,n] = size(Y); I = [];
        NanStats.NanMat = isnan(Y);
        
        NanStats.Feats = sum(NanStats.NanMat)*100/m;
        NanStats.Cases = sum(NanStats.NanMat,2)*100/n;
        
        NanStats.Cases25 = NanStats.Cases > 25 & NanStats.Cases <=50 ;
        NanStats.Feats25 = NanStats.Feats > 25 & NanStats.Feats <=50 ;
        NanStats.Cases25Perc = sum(NanStats.Cases25)*100/m;
        NanStats.Feats25Perc = sum(NanStats.Feats25)*100/n;
        
        NanStats.Cases50 = NanStats.Cases > 50 & NanStats.Cases <=99 ;
        NanStats.Feats50 = NanStats.Feats > 50 & NanStats.Feats <=99 ;
        NanStats.Cases50Perc = sum(NanStats.Cases50)*100/m;
        NanStats.Feats50Perc = sum(NanStats.Feats50)*100/n;
        
        NanStats.Cases100 = NanStats.Cases == 100;
        NanStats.Feats100 = NanStats.Feats == 100;
        NanStats.Cases100Perc = sum(NanStats.Cases100)*100/m;
        NanStats.Feats100Perc = sum(NanStats.Feats100)*100/n;
end
end

function Z = FillWithNan(Z, I)

        [~,n] = size(Z);
        tZ = nan(numel(I),n); 
        try
            tZ(~I,:) = Z;
        catch
            fprintf('error')
        end
        Z = tZ;

end