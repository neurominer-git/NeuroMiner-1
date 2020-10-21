function [ IN, OUT ] = SelectFeaturesAcrossCV1( IN, OUT, Param, minmaxfl, OptMode )
global VERBOSE

nc = size(OUT.F,3);
ns = size(IN.F,3);
nP = numel(Param.PFE.Perc);

funcomp = {@max, @min};

if check_fsizes(IN);
    cprintf('err','\nFeature masks have unequal dimensionalities due to feature preprocessing setup. Skipping Cross-CV1 feature operations.')
    return
end

qCrit = zeros(1,nP);
qIN=cell(1,nP); qOUT=cell(1,nP);
for q=1:nP
    
    qIN{q} = IN; qOUT{q} = OUT; qTrOpt = zeros(1,IN.nclass); qCVOpt = qTrOpt;
    
    for curclass=1:IN.nclass

        if VERBOSE, 
            if nP>1
                fprintf('\nAggregating feature masks for predictor #%g, threshold iteration: %g (T=%g)',curclass, q, Param.PFE.Perc(q)); 
            else
                fprintf('\nAggregating feature masks for predictor #%g',curclass); 
            end
        end
        kIndMat = [];

        for i=1:IN.nperms % loop through CV1 permutations

            for j=1:IN.nfolds  % loop through CV1 folds

                %% Aggregate feature subspace mask :
                if nc > 1 
                    % Every dichotomizer has its own feature subspace mask
                    kIndMat = [kIndMat qIN{q}.F{i,j,curclass}(:,qOUT{q}.kxVec{i,j,curclass}(qOUT{q}.F{i,j,curclass}))];
                elseif nc == 1 && ns ==1
                    kIndMat = [kIndMat qIN{q}.F{i,j}];
                else
                    % All dichotomizers share one feature subspace mask
                    % This mode is used for multi-group optimization
                    kIndMat = [kIndMat qIN{q}.F{i,j}(:,qOUT{q}.kxVec{i,j,curclass}(qOUT{q}.F{i,j}))];  
                end
            end
        end

        % Probabilistic Feature Extraction across CV1 partitions
        tkIndMat = ProbabilisticFea(kIndMat,Param.PFE, q);

        % Overwrite existing feature masks with probabilistic mask
        ll=1;
        for i=1:IN.nperms % loop through CV1 permutations

            for j=1:IN.nfolds  % loop through CV1 folds

                if size(tkIndMat,2) > 1, kIndMat_ll = tkIndMat(:,ll); else, kIndMat_ll = tkIndMat; end
                if ~any(kIndMat_ll),
                    if nP==1 || q==nP
                        error('\nYour feature selection procedure returned an empty feature space. Relax your selection parameters!'); 
                    else 
                        warning('\nYour feature selection procedure returned an empty feature space. Will check next threshold!'); 
                    end
                end
                if nc == 1 
                    if ns > 1
                        qIN{q}.F{i,j,curclass} = kIndMat_ll;
                    else
                        qIN{q}.F{i,j} = kIndMat_ll;
                    end
                else
                    qOUT{q}.F{i,j,curclass} = 1;
                    qOUT{q}.Weights{i,j,curclass} = 1;
                    qIN{q}.F{i,j,curclass} = kIndMat_ll;
                end    
                ll=ll+1;
            end
        end
    end

    % Retrain predictors using new feature masks
    [qIN{q}, qOUT{q}] = FoldPerm(qIN{q}, qOUT{q}, 'Retrain for PFC across CV1 partitions', OptMode, 0, 0, Param.SubSpaceStepping); 
    
    for curclass=1:IN.nclass
        qTrOpt(curclass) = nm_nanmean(nk_cellcat(qOUT{q}.tr(:,:,curclass),[],1))./nm_nansem(nk_cellcat(qOUT{q}.tr(:,:,curclass),[],1));
        qCVOpt(curclass) = nm_nanmean(nk_cellcat(qOUT{q}.ts(:,:,curclass),[],1))./nm_nansem(nk_cellcat(qOUT{q}.ts(:,:,curclass),[],1));
    end
    
    switch Param.datamode
        case 1
            qCrit(q) = mean(qTrOpt);
        case 2
            qCrit(q) = mean(qCVOpt);
        case 3
            qCrit(q) = mean([qTrOpt qCVOpt]);
    end
    
end
    
[opt, optq] = funcomp{minmaxfl}(qCrit);
if VERBOSE, fprintf('\tOptimum %g at T=%g', opt, Param.PFE.Perc(optq)); end
IN = qIN{optq}; OUT = qOUT{optq};
    
% Overwrite existing feature masks with probabilistic mask
for curclass=1:IN.nclass

    for i=1:IN.nperms % loop through CV1 permutations

        for j=1:IN.nfolds  % loop through CV1 folds

            OUT.TrHDperf(i,j,curclass) = OUT.tr{i,j,curclass};
            OUT.TrHTperf(i,j,curclass) = OUT.tr{i,j,curclass};       
            OUT.CVHDperf(i,j,curclass) = OUT.ts{i,j,curclass};
            OUT.CVHTperf(i,j,curclass) = OUT.ts{i,j,curclass}; 

        end
    end
end

% _________________________________________________________________________
% Helper function that checks whether dimensionalities of training data 
% partitions are identical or not
function [uneqfl, u_dims, max_u_dims] = check_fsizes(IN)

n_dims = zeros(IN.nperms * IN.nfolds,1);
cnt = 1;
for i=1:IN.nperms % loop through CV1 permutations
   for j=1:IN.nfolds  % loop through CV1 folds
       n_dims(cnt) = numel(IN.F{i,j});
       cnt=cnt+1;
   end
end

u_dims = unique(n_dims);

if numel(u_dims)>1, 
    uneqfl = true; 
    max_u_dims = max(u_dims);
else
    uneqfl = false;
    max_u_dims = u_dims;
end

