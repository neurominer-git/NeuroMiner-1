function [TransNodeMeanPerf, TransNodeSDPerf] = nk_ComputeTransNodePerformance(TransNodeArray, label, ix, jx, curclass)
global CV EVALFUNC MODEFL

[PermNum, FoldNum] = size(TransNodeArray);

Perf = zeros(PermNum, FoldNum);

for Perm = 1 : PermNum
    
    for Fold = 1 : FoldNum
        
        CVL = label(CV.TrainInd{ix,jx}(CV.cvin{ix,jx}.TestInd{Perm, Fold}));
        if strcmp(MODEFL,'classification')
            if ~exist('curclass','var') || isempty(curclass)
                ngroup = max(CVL);
                Prob = nk_ConvProbabilities(TransNodeArray{Perm,Fold}, ngroup);

                [mx, Pred] = max(Prob,[],2);
                Perf(Perm, Fold) = sum(Pred == CVL) / numel(CVL) * 100; 

            else

                Pred = sign(sum(TransNodeArray{Perm, Fold},2));
                if numel(CV.class{ix,jx}{curclass}.groups) > 1
                    indpos = CVL == CV.class{ix,jx}{curclass}.groups(1); indneg = CVL == CV.class{ix,jx}{curclass}.groups(2);
                    indall = CVL == CV.class{ix,jx}{curclass}.groups(1) | CVL == CV.class{ix,jx}{curclass}.groups(2);
                else
                    indpos = CVL == CV.class{ix,jx}{curclass}.groups;
                    indneg = ~indpos;
                    indall = true(size(indpos));
                end
                CVL(indpos) = 1; CVL(indneg) = -1; CVL(~indall) = [];
                Perf(Perm, Fold) = feval( EVALFUNC, CVL, Pred(indall,:));
            end
        else
            Pred = TransNodeArray{Perm, Fold};
            if size(Pred,2)>1, Pred = nm_nanmedian(Pred,2); end
            Perf(Perm, Fold) = feval( EVALFUNC, CVL, Pred);
        end
    end
    
end

TransNodeMeanPerf = mean(Perf(:));
TransNodeSDPerf = std(Perf(:));
