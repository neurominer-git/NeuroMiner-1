function [param, model] = nk_GetParam2(Y, label, Params, ModelOnly, FeatGroups)
% =========================================================================
% FORMAT [param, model] = nk_GetParam2(Y, label, Params, ModelOnly)
% =========================================================================
% Generic interface function to the algorithm-specific training modules
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2020
global TRAINFUNC SVM TIME CV CVPOS

% Remove cases which are completely NaN
[Y, label] = nk_ManageNanCases(Y, label);

timevec=[];
if ~isempty(TIME) && strcmp(SVM.prog,'WBLCOX')
    % Get training index
    TrInd = CV.TrainInd{CVPOS.CV2p,CVPOS.CV2f}(CV.cvin{CVPOS.CV2p,CVPOS.CV2f}.TrainInd{CVPOS.CV1p,CVPOS.CV1f});
    if CVPOS.fFull
        TrInd = [TrInd; CV.TrainInd{CVPOS.CV2p,CVPOS.CV2f}(CV.cvin{CVPOS.CV2p,CVPOS.CV2f}.TestInd{CVPOS.CV1p,CVPOS.CV1f})];
    end
    % Extract time vector
    timevec = TIME(TrInd);
    % Recode label into event [1] vs. no event [0]
    label(label==-1)=0;
end

% Run ADASYN if needed
if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
    [Y, label, timevec] = nk_PerfADASYN( Y, label, SVM.ADASYN, timevec); 
end

% Pass training matrix, labels, (and time vector) to used-defined training module
switch SVM.prog
    case 'SEQOPT'
        if  ~exist('FeatGroups','var') || isempty(FeatGroups)
            [param, model] = feval( TRAINFUNC, Y, label, [], ModelOnly, Params );
        else
            [param, model] = feval( TRAINFUNC, Y, label, FeatGroups, ModelOnly, Params );
        end
    case 'WBLCOX'
        [param, model] = feval( TRAINFUNC, Y, label, timevec, ModelOnly, Params );
    otherwise
        [param, model] = feval( TRAINFUNC, Y, label, ModelOnly, Params );
end

