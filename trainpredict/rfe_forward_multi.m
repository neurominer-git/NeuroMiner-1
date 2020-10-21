function [optparam, optind, optfound, optmodel] = ...
    rfe_forward_multi(Y, mY, label, labelB, labelM, Ynew, labelnew, labelnewM, Ps, FullFeat, FullParam, ngroups, ActStr)
%==========================================================================
% [optparam, optind, optfound, optmodel] = ...
%   rfe_forward(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr)
%==========================================================================
% Greedy forward search wrapper algorithm
%--------------------------------------------------------------------------
% INPUT:
%  Y             :   CV1 Training data
%  label         :   CV1 Training labels
%  Ynew          :   CV1 Test data
%  labelnew      :   CV1 Test labels
%  Ps            :   Parameters for learning algorithm
%  FullFeat      :   Feature mask prior to wrapper-based feat extraction
%  FullParam     :   Model performance associated with FullFeat
%  ActStr        :   'Tr'   use CV1 training data performance, 
%                    'Ts'   use CV1 test data performance,
%                    'TrCV' use both CV1 training and CV1 test data to
%                    compute performance
%
% OUTPUT:
%  optparam      :   Performance of wrapper-optimized feature set
%  optind        :   Indices to selected features
%  optfound      :   Flag whether an optimal model has been identified
%  optmodel      :   Optimized model
%==========================================================================
%(c) Nikolaos Koutsouleris, 05/2017
global VERBOSE TRAINFUNC 

r = rfe_algo_settings_multi(Y, mY, label, labelB, labelM, Ynew, labelnew, labelnewM, Ps, FullFeat, FullParam, ngroups, ActStr);
nclass = numel(Y);
optfound = 0; optparam = r.optparam; r.InitializeOrder = true;

if VERBOSE
    fprintf('\n-----------------------------')
    fprintf('\nGREEDY FORWARD FEATURE SEARCH')
    fprintf('\n-----------------------------')
    fprintf('\nOptimization data mode: %s', ActStr)
    fprintf('\nParameter evaluation: %s (%s)', r.evaldir, r.optfunc)
    if r.FeatStepPerc
        for curclass=1:nclass
            fprintf('\nStepping: %g%% of %g features in model #%g per wrapper cycle.',r.lperc(curclass), r.kFea(curclass), curclass)
        end
    else
        fprintf('\nStepping: Top feature in wrapper cycle.')
    end
    if r.FeatRandPerc
        fprintf('\nRandom feature selection: %g%% of top-ranked features in block',r.FeatRandPerc)
    end
end

 %% Start Wrapper: FORWARD FEATURE SELECTION
S = cell(1,nclass); Sind = cell(1,nclass); k = zeros(1,nclass); lstep = zeros(1,nclass);
for curclass=1:nclass
    k(curclass) = r.kFea(curclass); 
    if r.PreSort
        Sind{curclass} = r.PreOrder{curclass}; 
    else
        Sind{curclass} = 1:r.kFea(curclass); 
    end
    lstep(curclass) = 1; if r.lperc(curclass), lstep(curclass) = ceil((numel(Sind{curclass})/100)*r.lperc(curclass)); end
end
maxK = max(k); minK = min(r.MinNum); 
Opt = struct('S',[],'Param',[],'ParamTs',[]);
cnt=1;     

while maxK > minK

    maxNS = max(cellfun(@numel,Sind));
    if ~maxNS, break; end
    val = zeros(maxNS,1); 

    if VERBOSE, 
        for curclass=1:nclass,fprintf('\n\tFeature pool size of model #%g: %4.0f out of %4.0f, block size: %4.0f feature(s) ', ...
                curclass, numel(S{curclass}), numel(Sind{curclass}), lstep(curclass)); end
    end
    
    tEnd = zeros(maxNS, nclass);
    ds = zeros(size(r.T{1},1), nclass, maxNS );
    ts = zeros(size(r.T{1},1),nclass, maxNS );
    while maxNS > 0

        for curclass=1:nclass
            NSind = numel(Sind{curclass});
            if maxNS > NSind, 
                tEnd(maxNS, curclass) = Sind{curclass}(NSind); 
            else, 
                tEnd(maxNS, curclass) = Sind{curclass}(maxNS);
            end 
            kS = [S{curclass} tEnd(maxNS,curclass)]; 
            tY = r.Y{curclass}(:,kS); T = r.T{curclass}(:,kS);
            [~, model] = feval(TRAINFUNC, tY, label{curclass}, 1, Ps{curclass});
            [~, ts(:, curclass, maxNS), ds(:, curclass, maxNS)] = nk_GetTestPerf(T, r.L{curclass}, [], model, tY);
        end
        val(maxNS) = nk_MultiEnsPerf(ds(:,:, maxNS), ts(:,:,maxNS), r.Lm, 1:nclass, r.ngroups);
        maxNS = maxNS - 1;
    end

     % Sort features
    if r.perm
        [valperm, permind] = rfe_multi_perm_featcombs(ds, ts, r.Lm, nclass,r.nperms, r.ngroups);
        [~, ind] = sort(valperm(:), r.optfunc);
    else
        [~, ind] = sort(val,r.optfunc);
    end
    ds = zeros(size(r.T{1},1),nclass);
    ts = zeros(size(r.T{1},1),nclass);
    krem = 1:min(lstep);rind = krem;

    for curclass=1:nclass
        if r.perm
            IndAdd = tEnd(permind(curclass,ind(krem)),curclass)';
        else
            IndAdd = tEnd(ind(krem),curclass)';
        end
        kS = [S{curclass} IndAdd];
        tY = r.Y{curclass}(:,kS); T = r.T{curclass}(:,kS);
        [~, model] = feval(TRAINFUNC, tY, label{curclass}, 1, Ps{curclass});
        [~, ds(:,curclass), ts(:,curclass)] = nk_GetTestPerf(T, r.L{curclass}, [], model, tY); 
    end
    param = nk_MultiEnsPerf(ds, ts, r.Lm, 1:nclass, r.ngroups);      

    % Add feature to feature space only if current performance is better
    % then previous space
    % Finally trained model will have only a performance == optparam
    % if the feature stepping is 1
    if feval(r.evaldir, param, optparam) 
        % If user activated random feature selection then we have
        % to randomly choose x% of the features in the block
        if r.FeatRandPerc
            rstep = ceil(min(lstep)/100)*r.FeatRandPerc;
            rind = randperm(min(lstep),rstep);
        % ... otherwise select entire block
        else
            rind = krem;
        end
        for curclass=1:nclass
            if r.perm
                IndAdd = tEnd(permind(curclass,ind(rind)),curclass)';
            else
                IndAdd = tEnd(ind(krem),curclass)';
            end
            S{curclass} = [S{curclass} IndAdd]; 
            Opt.S{curclass, cnt} = S{curclass};
        end
        optparam = param; optfound = 1; 
        if VERBOSE, fprintf('=> NEW optimum: # Features: %4.0f ==> %s = %g', mean(cellfun(@numel,S)), ActStr, optparam); end
        Opt.Param = [Opt.Param optparam];
        %Opt.ParamTs = [Opt.ParamTs param_ts];
        cnt=cnt+1;
    end

    for curclass=1:nclass
        if r.perm
           IndAdd = tEnd(permind(curclass,ind(rind)),curclass)';
        else
           IndAdd = tEnd(ind(krem),curclass)';
        end
        remInd = ismember(Sind{curclass},IndAdd);
        Sind{curclass}(remInd) = [];
        % Recompute lstep according to current feature pool
        if r.lperc(curclass), lstep(curclass) = ceil((numel(Sind{curclass})/100)*r.lperc(curclass)); end
    end
    maxK = maxK - numel(rind);

end
        
for curclass=1:nclass
    if r.KneePoint,
        kneepoint = knee_pt(Opt.Param,[],true);
        if VERBOSE,
            if isnan(kneepoint)
                cprintf('err','\nNot enough data points to compute kneepoint. Selecting final feature mask.');
            else
                fprintf('\nSelected kneepoint of optimization curve at wrapper cycle #%g => %s = %g', kneepoint, ActStr, Opt.Param(kneepoint));
            end
        end
        if isnan(kneepoint), kneepoint = numel(Opt.S); end
        optind = r.FullInd{curclass}(Opt.S{curclass}{kneepoint});

    elseif isnan(param) && ~optfound 
        % Optimization returned non-finite performance at given parameter
        % combination. Return original feature space
        optind{curclass} = r.FullInd{curclass};
        fprintf('\n');warning('Greedy forward search did not return any feature mask for given parameter setting. Return original feature space.')
    else
        optind{curclass} = r.FullInd{curclass}(S{curclass});
    end
    if VERBOSE, fprintf('\nDone. '); end
    [~,optmodel{curclass}] = feval(TRAINFUNC, Y{curclass}(:,optind{curclass}), label{curclass}, 1, Ps{curclass}); 
end
optfound = 1;

