function [R, optmodel] = nk_MLOptimizer_WrapperMulti(Y, mY, L, bL, MultiL, mYnew, Lnew, MultiLnew, ngroups, Ps, FullParam, SubFeat)

global RFE SVM

if nargin < 12, 
    for curclass=1:numel(Y)
        SubFeat{curclass} = true(1,size(Y{curclass},2));
    end
end
ActStr = {'Tr', 'CV', 'TrCV'};

% Remove cases which are completely NaN
for i=1:numel(Y)
    [Y{i}, L{i}] = nk_ManageNanCases(Y{i}, L{i});
    
    % Run ADASYN if needed
    if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
        [Y{i}, L{i}] = nk_PerfADASYN(Y{i}, L{i}, SVM.ADASYN);
    end
    
    if i==1
        [mY{i}, MultiL, I] = nk_ManageNanCases(mY{i}, MultiL); 
        [mYnew{i}, Lnew{i}, Inew] = nk_ManageNanCases(mYnew{i}, Lnew{i}); 
        MultiLnew(Inew)=[]; bL{i}(I)=[];
    else
        mY{i}(I,:)=[]; bL{i}(I)=[];
        mYnew{i}(Inew,:)=[]; Lnew{i}(Inew)=[]; 
    end
    
end

switch RFE.Wrapper.type
    %% GREEDY FORWARD/BACKWARD FEATURE SEARCH
    case 1 
        funs = { @rfe_forward_multi,  @rfe_backward_multi };
        [optparam, optind, optfound, optmodel] = funs{RFE.Wrapper.GreedySearch.Direction}(Y, mY, L, bL, MultiL, mYnew, Lnew, MultiLnew, Ps, SubFeat, FullParam, ngroups, ActStr{RFE.Wrapper.datamode});
    
    %% SIMULATED ANNEALING
    % [ this does not work in multi-class mode yet ]
    case 2
        [optparam, optind, optfound, optmodel] = nk_SimAnneal(Y, label, Ynew, labelnew, Ps, SubFeat, FullParam, ActStr{RFE.Wrapper.datamode});
end
% Transfer params to output structure
R.found                   = optfound;
R.FeatureIndex            = optind;
if nargin == 12
    for curclass=1:nclass
        R.SubFeatIndex{curclass}                   = false(size(SubFeat{curclass}));
        R.SubFeatIndex{curclass}(optind{curclass}) = SubFeat(optind{curclass}); 
    end
end
R.OptimParam              = optparam;


