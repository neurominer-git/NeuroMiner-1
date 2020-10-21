function [MultiGroupSelection, MultiBinBind] = nk_MultiBinBind(GD, label, f, d, curlabel)
global CV RFE

%%%% PREPARATIONS %%%%
ngroups                 = numel(unique(label(~isnan(label))));
nclass                  = length(CV.class{1,1});
[CV1perms, CV1folds]    = size(CV.cvin{f,d}.TrainInd);
TsXdims                 = size(CV.TestInd{f,d},1);
nc                      = size(GD.BinaryGridSelection{1}{curlabel}.bestfeats,3);
Ydims = zeros(nclass,1); Npos = cell(nclass,1); 
MultiGroupSelection.SelNodes = GD.BinaryGridSelection{1}{curlabel}.SelNodes;
for curclass=1:nclass
    if nc > 1
        Ydims(curclass) = nk_GetCV2EnsembleDims(GD.BinaryGridSelection{curclass}{curlabel}.bestfeats(:,:,curclass));
    else
        Ydims(curclass) = nk_GetCV2EnsembleDims(GD.BinaryGridSelection{curclass}{curlabel}.bestfeats);
    end
    Npos{curclass} = GD.BinaryGridSelection{curclass}{curlabel}.Npos;
    MultiGroupSelection.bestP{curclass} = GD.BinaryGridSelection{curclass}{curlabel}.bestP;
    IndSelNodes = find(GD.BinaryGridSelection{curclass}{curlabel}.SelNodes);
    MultiGroupSelection.SelNodes(IndSelNodes) = true; 
    MultiGroupSelection.bestfeats{curclass} = cell(CV1perms, CV1folds,numel(Npos{curclass}));
end
mTs                     = zeros(TsXdims,sum(Ydims));
Classes                 = zeros(1,sum(Ydims));
MultiBinBind.TrPerf     = zeros(CV1perms,CV1folds);
MultiBinBind.CVPerf     = zeros(CV1perms,CV1folds);
MultiBinBind.Ts1Perf    = zeros(CV1perms,CV1folds);
mcolend                 = 0;
%%%% CV1 loop
for k=1:CV1perms

    for l=1:CV1folds
        
        mTrLabel = label(CV.TrainInd{f,d}(CV.cvin{f,d}.TrainInd{k,l}));
        mCVLabel = label(CV.TrainInd{f,d}(CV.cvin{f,d}.TestInd{k,l}));
        mTsLabel = label(CV.TestInd{f,d});
        %TrXdims     = size(CV.cvin{f,d}.TrainInd{k,l},1);
        %CVXdims     = size(CV.cvin{f,d}.TestInd{k,l},1);
        mTr         = [];
        mCV         = [];

        %%%% DICHOTOMIZATION LOOP 
        for curclass = 1:nclass
            
            MultiGroupSelection.bestcomplexity(curclass) = mean(GD.C(Npos{curclass}));
            
            for curnode = 1:numel(Npos{curclass})
            
                Fkl     = GD.FEAT{Npos{curclass}(curnode)};
                if size(Fkl,3) > 1, Fkl = Fkl{k,l,curclass}; else Fkl = Fkl{k,l}; end
                MultiGroupSelection.bestfeats{curclass}{k,l,curnode} = Fkl;
                
                ul      = size(Fkl,2);

                % Get decision values at binary classifiers best positions
                TrD    = GD.DT{Npos{curclass}(curnode)}{k,l,curclass};
                CVD    = GD.DV{Npos{curclass}(curnode)}{k,l,curclass}; 
                TsD    = GD.DS{Npos{curclass}(curnode)}{k,l,curclass};

                % Compute column pointers for multi-group CV2 array construction
                mcolstart = mcolend + 1; mcolend = mcolstart + (ul-1);
                if curclass == 1 && curnode ==1 , mcolX = mcolstart; end;

                % Enter predictions for multi-group classification into arrays
                mTr = [mTr TrD];
                mCV = [mCV CVD];
                mTs(:,mcolstart:mcolend) = TsD;
                Classes(1,mcolstart:mcolend) = curclass;
            end
        end
        
        % Compute multi-group performance on CV1 training data for CV1
        % partition [k,l]
        [MultiBinBind.TrPerf(k,l), MultiBinBind.TrPred{k,l}] = ...
            nk_MultiEnsPerf(mTr, sign(mTr), mTrLabel, Classes(:,mcolX:mcolend), ngroups);
                        
        % Compute multi-group performance on CV1 test data for CV1
        % partition [k,l]
        [MultiBinBind.CVPerf(k,l), MultiBinBind.CVPred{k,l}] = ...
            nk_MultiEnsPerf(mCV, sign(mCV), mCVLabel, Classes(:,mcolX:mcolend),ngroups);
                        
        % Compute multi-group performance on CV2 test data for CV1
        % partition [k,l]                        
        [MultiBinBind.Ts1Perf(k,l), MultiBinBind.Ts1Pred{k,l}] = ...
            nk_MultiEnsPerf(mTs(:,mcolX:mcolend), ...
                            sign(mTs(:,mcolX:mcolend)), ...
                            mTsLabel, ...
                            Classes(:,mcolX:mcolend), Groups);
    end
end

MultiBinBind.Mean_TrPerf    = mean(MultiBinBind.TrPerf(:));
MultiBinBind.SD_TrPerf      = std(MultiBinBind.TrPerf(:));
MultiBinBind.Mean_CVPerf    = mean(MultiBinBind.CVPerf(:));
MultiBinBind.SD_CVPerf      = std(MultiBinBind.CVPerf(:));
MultiBinBind.Mean_Ts1Perf   = mean(MultiBinBind.Ts1Perf(:));
MultiBinBind.SD_Ts1Perf     = std(MultiBinBind.Ts1Perf(:));

% Compute multi-group performance on CV2 test data by using entire CV1
% prediction data on CV2 instances (ensemble of ensembles decision)
[MultiBinBind.Ts2Perf, MultiBinBind.Ts2Pred] = nk_MultiEnsPerf(mTs, sign(mTs), mTsLabel, Classes, ngroups);

MultiGroupSelection.bestacc = MultiBinBind.Mean_CVPerf;
switch RFE.CV2Class.type
    case 1
        MultiGroupSelection.besttestparam  = MultiBinBind.Mean_Ts1Perf;
    case 2
        MultiGroupSelection.besttestparam  = MultiBinBind.Ts2Perf;
end

MultiGroupSelection.bestcomplexity    = mean(MultiGroupSelection.bestcomplexity);
MultiGroupSelection.binbind           = true;
MultiGroupSelection.bestpred          = MultiBinBind.Ts2Pred;
MultiGroupSelection.bestCV2pred       = MultiBinBind.Ts1Pred;

% Print some info
fprintf('\nMulti-group performance at dichotomizers'' optima:\nCV1 = %1.2f, CV2 = %1.2f, Complexity = %1.2f', ...
    MultiGroupSelection.bestacc, MultiGroupSelection.besttestparam, MultiGroupSelection.bestcomplexity)
for curclass=1:nclass
    for curnode = 1:numel(Npos{curclass})
        if numel(Npos{curclass})>1
            fprintf('\nParameter combination at node #%g (%s):\t', curnode, CV.class{1,1}{curclass}.groupdesc)
            fprintf('%5.2f', GD.BinaryGridSelection{curclass}{curlabel}.bestP(curnode,:))
        else
            fprintf('\nParameter combination (%s):\t', CV.class{1,1}{curclass}.groupdesc)
            fprintf('%5.2f', GD.BinaryGridSelection{curclass}{curlabel}.bestP)
        end
    end
end
