% =========================================================================
% function [P, Perf, d, crit] = nk_MultiDecideHierarchOneVsOne(H, Y,
%                                                   Classes, maxfunc)
% =========================================================================
% 
% INPUTS:
% -------
% H :                   Ensemble consisting of binary one-vs-one dichotomizer
%                       prediction for given observations
% Y :                   Multi-group label vector
% Classes :             Dichotomization vector with 1 -> number of dichotomizers in
%                       ensemble
% maxfunc :             Decision function, with 1 = sum, 2 = mean, 3 = product, 
%                       4 = majority vote
%
% OUTPUTS:
% --------
% P :                   Multi-group predictions
% Perf :                Multi-group prediction performance
% d :                   decision matrix
% crit :                selected first-rank (max) decision score from d
%
% COMMENTS:
% ---------
% This function performs Hierarchical One-Vs-One-Max-Wins multi-group 
% classification. Initially, the first group is tested against
% the rest, then the second against the rest, until the final group is
% defined by the last binary dichotomizer.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2011

function [P, Perf, d, crit] = nk_MultiDecideHierarchOneVsOne(H, Y, Classes, maxfunc)

global SVM

if ~isempty(SVM) && SVM.GridParam == 14
        multimode = 1;
    else
        multimode = 0;
end

if iscell(H)
    [iy,jy,nclass]  = size(H);
    P               = cell(iy,jy);
    crit            = cell(iy,jy);
    Perf            = zeros(iy,jy);
else
    iy = 1; jy=1;
    nclass          = max(Classes);
end

for k=1:iy
    for l=1:jy
        if iscell(H)
            nsubj       = size(H{k,l},1);
            M           = nk_OneVsOne(Classes{k,l}, nclass);
            ngroups     = size(M,1);
            P{k,l}      = zeros(nsubj,1);
        else
            nsubj       = size(H,1);
            M           = nk_OneVsOne(Classes, nclass);
            ngroups     = size(M,1);
        end
        
        td              = zeros(nsubj,1);
        
        % Perform Hierarchival One-vs-Rest
        % Step 1) 1st group against the rest
        if iscell(H),
            tH= [];
            for curclass=1:nclass
                tH = [tH H{k,l,curclass}]; 
            end
        else
            tH = H; 
        end

        for j=1:ngroups-1
            
            % Get codeword for current group
            mOne          = M(j,:) > 0;
            
            % Get index to remaining data points without label
            indSubj       = find(td==0);
            tHJ           = tH(indSubj,mOne); 
            
            switch maxfunc
                case 1
                    gOne = sum(tHJ,2); 
                case 2
                    gOne = mean(tHJ,2); 
                case 3
                    gOne = prod(tHJ,2); 
                case 4
                    gOne = sum(sign(tHJ),2); 
            end
            
            % Decide
            td(indSubj(gOne>0)) = j;
        end
        
        td(td==0) = ngroups;
        
        if iscell(H)
            P{k,l} = td;
            if isempty(Y), Y = ones(numel(P{k,l}),1); end                       
            Perf(k,l) = nk_MultiPerfQuant(Y, P{k,l}, multimode);
        else
            P = td;
            if isempty(Y), Y = ones(numel(P),1); end
            Perf(k,l) = nk_MultiPerfQuant(Y, P, multimode);
        end
        
    end
end

return
