function [CatArray, ClassArray] = nk_CatNodes(arr, curclass)


NodeNum = numel(arr);

[PermNum, FoldNum, nclass] = size(arr{1});

CatArray = cell(PermNum, FoldNum);
ClassArray = cell(PermNum, FoldNum);

for Perm = 1:PermNum
        
    for Fold = 1:FoldNum

        CatArray{Perm, Fold} = [];
        ClassArray{Perm, Fold} = [];
        
        for Node = 1 : NodeNum
            
            if ~exist('curclass','var') || isempty(curclass)
                try
                    CatArray{Perm, Fold} = [ CatArray{Perm, Fold} arr{Node}{Perm, Fold} ];
                catch
                    CatArray{Perm, Fold} = [ CatArray{Perm, Fold} arr{Node} ];
                end
            else
                CatArray{Perm, Fold} = [ CatArray{Perm, Fold} arr{Node}{Perm, Fold, curclass} ];
            end
        end

    end

end
               