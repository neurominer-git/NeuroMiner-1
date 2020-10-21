%==========================================================================
%FORMAT [F, Weights, Mx] = EvalSubSpaces(Perf, strout, ...             
%                                   SubSpaceStrategy, EnsStrat, ...     
%                                   Crit, minmaxfl, AddPerf, Weighting) 
%==========================================================================
%'EVALUATE SUBSPACES' MODULE                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) Nikolaos Koutsouleris, 08/2018

function [F, Weights, Mx] = EvalSubSpaces(Perf, strout, SubSpaceStrategy, Crit, minmaxfl, AddPerf, Weighting, Pred)

global VERBOSE

[nperms, nfolds, nclass] = size(Perf);
F       = cell(nperms, nfolds, nclass);
Mx      = zeros(nperms, nfolds, nclass);
Weights = cell(nperms, nfolds, nclass);

% Check whether each fold contains equal number of subspaces
sz = zeros(numel(Perf),1);
for i=1:numel(Perf)
    sz(i) = size(Perf{i},1);
end
if numel(unique(sz)) > 1,
    avgflag = false;
    if VERBOSE; fprintf('\nUnequal number of subspaces across CV1 data partitions.'); end
else
    avgflag = true;
end

if avgflag
    if VERBOSE; fprintf('\nDetermine subspace optimum across folds & permutations.'); end
    PerfAvg = cell(nclass,1);
    for curclass=1:nclass
        PerfAvg{curclass} = mean(nk_cellcat(Perf(:,:,curclass),[],2),2);
    end
else
    if VERBOSE; fprintf('\nDetermine subspace optimum in each fold separately.'); end
end

for i=1:nperms  
    
    for j=1:nfolds
        
        for curclass=1:nclass
            
            if avgflag
                Px = PerfAvg{curclass};
            else
                Px = Perf{i,j,curclass};
            end
            
            kFea = size(Perf{i,j,curclass},1);
            switch SubSpaceStrategy 
                case 0
                    F{i,j,curclass} = true;
                case 4
                    F{i,j,curclass} = true(kFea,1);
                otherwise
                    if VERBOSE, fprintf('\n%s => CV1 [perm %g, fold %g]:\t',strout,i,j); end
                    F{i,j,curclass} = false(kFea,1);
                    switch minmaxfl

                        case 1 % Maximize
                            % Maximize ranking criterion across feature subspaces
                            % Additional performance array is provided in order to
                            % avoid overfitting
                            if ~isempty(AddPerf) 
                                inx = Px <= AddPerf{i,j,curclass};
                                if ~any(inx)
                                    [Mx(i,j,curclass),indm] = min(Px);
                                else
                                    [Mx(i,j,curclass),indm] = max(Px(inx));
                                end                  
                            else
                                [Mx(i,j,curclass),indm] = max(Px);
                            end

                            switch SubSpaceStrategy
                                case 1 % "the winner takes it all"              
                                    F{i,j,curclass}(indm) = true;
                                    if VERBOSE, fprintf('max=%1.4f',Mx(i,j,curclass)); end
                                case 2 % within range defined by user from max
                                    F{i,j,curclass}( Px >= (Mx(i,j,curclass) - Crit) ) = true;
                                    if VERBOSE, fprintf('max=%1.4f, min=%1.4f',Mx(i,j,curclass),(Mx(i,j,curclass) - Crit)); end
                                case 3 % percentile defined by user
                                    pr = percentile( Px, Crit ) ;
                                    F{i,j,curclass}( Px >= pr ) = true;
                                    if VERBOSE, fprintf('max=%1.4f, min=%1.4f',Mx(i,j,curclass), pr); end;
        %                       case 4 % adaptive method (ROC for classification/MSE for regression problems)
        %                             % Sort feature subspaces
        %                             [sPx,sInd] = sort(Px,'descend');
        %                             PerfEns = zeros(kFea,1);
        %                             for z=1:kFea
        %                                 % Build ensemble
        %                                 EnsDat = D{i,j}(:,sInd(1:z));
        %                                 % Evaluate ensemble
        %                                 switch MODEFL
        %                                     case 'classification'
        %                                         JointHypo = sign(mean(EnsDat,2));
        %                                         PerfEns(z) = get_test_param(SVM.GridParam,CV)
        %                                     case 'regression'
        %                                         JointHypo = mean(EnsDat,2)
        %                                 end
        %                             end
                            end

                        case 2 % Minimize
                            % Minimize ranking criterion across feature subspaces
                            [Mx(i,j,curclass),indm] = min(Px);
                            switch SubSpaceStrategy
                                case 1 % "the winner takes it all"              
                                    F{i,j,curclass}(indm) = true;
                                     if VERBOSE, fprintf('max=%1.4f',Mx(i,j,curclass)); end;
                                case 2 % within range defined by user from max
                                    F{i,j,curclass}( Px <= (Mx(i,j,curclass) + Crit) ) = true;
                                    if VERBOSE, fprintf('max=%1.4f, min=%1.4f',Mx(i,j,curclass),(Mx(i,j,curclass) - Crit)); end;
                                case 3 % percentile defined by user
                                    pr = percentile(Px,Crit);
                                    F{i,j,curclass}( Px <= pr ) = true;
                                    if VERBOSE, fprintf('max=%1.4f, min=%1.4f',Mx(i,j,curclass), pr); end;
                            end
                    end
            end
            
            indm = F{i,j,curclass};
            %Weights{i,j,curclass} = zeros(size(indm));
            
            switch Weighting
                case 0 % equal weights
                    Weights{i,j,curclass} = ones(size(indm));
                case 1 % compute linear weighting
                    Weights{i,j,curclass} = Px(indm)/mean(Px(indm));
                case 2 % compute logarithmic weighting    
                    Weights{i,j,curclass} = log(Px(indm))/mean(log(Px(indm)));
            end
        end
    end
end

