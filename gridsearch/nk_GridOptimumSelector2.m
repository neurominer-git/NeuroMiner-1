function [GridSelectorResults, ipos] = ...
    nk_GridOptimumSelector2(TR, TS, C, MD, FEAT, Weights, CV2Pred, CV1Pred, P, Pdesc, combcell, act, perc)

global SAV 

if ~exist('perc','var'), perc = []; end

%%%%%%%%%%%%%%%%%%%%%%%% MAX SELECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
C(isnan(C)) = 0;
[ipos, ind0] = nk_FindGridOpt2(TR, TS, C, act, perc);
nipos = numel(ipos); 
%%%%%%%%%%%%%%%%%%% SELECT PARAMS AT OPTIMUM %%%%%%%%%%%%%%%%%%%%
if iscell(P) 
    sP = size(P);
    if ~combcell
        for curclass=1:numel(P)
           GridSelectorResults.bestP{curclass}= P{curclass}(ipos,:);
        end
    else
        GridSelectorResults.bestP = P(ipos,:);
    end
else
    GridSelectorResults.bestP = P(ipos,:);
end

nP = 1;
if iscell(P) && ~combcell, 
    nP = sP(2); 
elseif ~iscell(P) && (numel(P)==1 && isnan(P))
    P = repmat({nan},nP,1); Pdesc = P;
end

for curclass=1:nP
    if (iscell(P) && ~combcell) || iscell(Pdesc{1}) 
        tPdesc = Pdesc{curclass}; tP = P{curclass}; 
    else
        tPdesc = Pdesc; tP = P; 
    end
    if ~isempty(tPdesc)
        headstr = sprintf('\nNo. ');
        for n=1:numel(tPdesc)
            if iscell(tPdesc)
                headstrn = sprintf('%15s',tPdesc{n});
            else
                headstrn = sprintf('%15s',tPdesc(n));
            end
            headlen(n) = numel(headstrn);
            headstr = sprintf('%s|%s',headstr,headstrn);
        end
        headstr = sprintf('%s |%12s |%13s\n',headstr, 'CV1', 'CV2');
        cprintf('*blue','%s',headstr);
        cprintf('*blue',repmat('=',1,length(headstr)));
        for i=1:nipos
            paramstr = [];
            for n = 1:size(tP(ipos(i),:),2)
               if combcell
                   if isnumeric(tP{ipos(i),n})
                       tPn = num2str(tP{ipos(i),n});
                   else
                       tPn = tP{ipos(i),n};
                   end
               else
                   tPn = num2str(tP(ipos(i),n));  
               end
               paramstr = [ paramstr sprintf(['%' num2str(headlen(n)) 's|'], tPn ) ];
            end
            fprintf('\n[%4g]%s%13.2f|%13.2f', i, paramstr, TR(ipos(i)), TS(ipos(i)))
        end 
        fprintf('\n');
        cprintf('*blue',repmat('=',1,length(headstr)));
    end
 end


GridSelectorResults.Npos = ipos;
GridSelectorResults.SelNodes = ind0;

if SAV.savemodel, GridSelectorResults.bestmodel = MD(ipos); end
GridSelectorResults.bestfeats  = FEAT(ipos);
GridSelectorResults.bestweights = Weights(ipos);
GridSelectorResults.bestacc = TR(ipos);
GridSelectorResults.bestpred = CV2Pred(ipos);
GridSelectorResults.bestCV1TsPred = CV1Pred(ipos);
GridSelectorResults.besttestparam = TS(ipos);
GridSelectorResults.besterr = TS(ipos) - TR(ipos);

fprintf('\nMean CV1 performance: %1.2f, Mean CV2 performance: %g across %g nodes', ...
    mean(GridSelectorResults.bestacc), mean(GridSelectorResults.besttestparam), ...
    numel(GridSelectorResults.Npos))

if ~isempty(C),
    GridSelectorResults.bestcomplexity = C(ipos);
end
GridSelectorResults.Nodes = length(ipos);

end