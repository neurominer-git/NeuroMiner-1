function [optparam, optind, optfound, optmodel] = ...
    rfe_backward_multi(Y, mY, label, labelB, labelM, Ynew, labelnew, labelnewM, Ps, FullFeat, FullParam, ngroups, ActStr)
%=================================================================================================
% [optparam, optind, optfound, optmodel] = ...
%   rfe_backward(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr)
%=================================================================================================
% Greedy recursive feature elimination algorithm
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
global VERBOSE BATCH TRAINFUNC 

r = rfe_algo_settings_multi(Y, mY, label, labelB, labelM, Ynew, labelnew, labelnewM, Ps, FullFeat, FullParam, ngroups, ActStr);
nclass=numel(Y);

%% Feature block size settings
optfound = 0; optparam = r.FullParamMulti; optind = r.FullInd;
S = cell(1,nclass); k = zeros(1,nclass); lstep = zeros(1,nclass);
for curclass=1:nclass
    k(curclass) = r.kFea(curclass); 
    if r.PreSort
       S{curclass} = r.PreOrder{curclass};
    else
       S{curclass} = 1:r.kFea(curclass); 
    end
    lstep(curclass) = ceil((numel(S{curclass})/100)*r.lperc(curclass));
end
optmodel = r.FullModel;

if max(k) <= r.MinNum, return; end

Opt = struct('S',[],'Param',[],'ParamTs',[]);

%% Start Wrapper: RECURSIVE FEATURE ELIMINATION
if VERBOSE
    fprintf('\n-----------------------------------')
    fprintf('\nGREEDY BACKWARD FEATURE ELIMINATION')
    fprintf('\n-----------------------------------')
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

maxK = max(k); minK = min(r.MinNum);
cnt=1;
    
while maxK > minK

    maxNS = max(cellfun(@numel,S));
    if ~maxNS, break; end
    val = zeros(maxNS,1); 

    if VERBOSE, 
        for curclass=1:nclass,fprintf('\n\tFeature pool size of model #%g: %4.0f, block size: %4.0f feature(s) ', ...
                curclass, numel(S{curclass}), lstep(curclass)); end
    end

    tEnd = zeros(maxNS, nclass);
    ds = zeros(size(r.T{1},1), nclass, maxNS );
    ts = zeros(size(r.T{1},1), nclass, maxNS );

    while maxNS > 0

        for curclass=1:nclass

            NS = numel(S{curclass});
            if maxNS > NS, 
                tEnd(maxNS, curclass) = NS; 
            else, 
                tEnd(maxNS, curclass) = maxNS;
            end 
            lind = 1:numel(S{curclass}); lind(tEnd(maxNS,curclass))=[]; kS = S{curclass}(lind);
            tY = r.Y{curclass}(:,kS); T = r.T{curclass}(:,kS);
            [~, model] = feval(TRAINFUNC, tY, label{curclass}, 1, Ps{curclass});
            [~, ds(:,curclass, maxNS), ts(:,curclass, maxNS)] = nk_GetTestPerf(T, r.L{curclass}, [], model, tY);
        end
        val(maxNS) = nk_MultiEnsPerf(ds(:,:, maxNS), ts(:,:,maxNS), r.Lm, 1:nclass, r.ngroups);
        maxNS = maxNS - 1;
    end

    % Eliminate feature, that maximally degrades CV1-Test performance
    if r.perm
        [valperm, permind] = rfe_multi_perm_featcombs(ds, ts, r.Lm, nclass,r.nperms, r.ngroups);
        [~, ind] = sort(valperm(:), r.optfunc);
    else
        [~, ind] = sort(val,r.optfunc);
    end
    
    krem = 1:min(lstep);

    % If user activated random feature selection then we have
    % to randomly choose x% of the features in the block
    if r.FeatRandPerc
        rstep = ceil(min(lstep)/100)*r.FeatRandPerc;
        rind = randperm(min(lstep),rstep);
    % ... otherwise select entire block
    else
        rind = krem;
    end
    
    ds = zeros(size(r.T{1},1),nclass);
    ts = zeros(size(r.T{1},1),nclass);
    
    for curclass=1:nclass
        if r.perm
            remInd = ismember(S{curclass},tEnd(permind(curclass,ind(rind)),curclass));
        else
            remInd = ismember(S{curclass},tEnd(ind(rind),curclass));
        end
        S{curclass}(remInd) = [];
        lstep(curclass) = ceil((max(cellfun(@numel,S))/100)*min(r.lperc));
        
        tY = r.Y{curclass}(:,S{curclass}); T = r.T{curclass}(:,S{curclass});
        [~, model] = feval(TRAINFUNC, tY, label{curclass}, 1, Ps{curclass});
        [~, ds(:,curclass), ts(:,curclass)] = nk_GetTestPerf(T, r.L{curclass}, [], model, tY); 
    end
    param = nk_MultiEnsPerf(ds, ts, r.Lm, 1:nclass, r.ngroups);   

    if feval(r.evaldir, param, optparam)
        optfound = 1; optparam = param; 
        for curclass=1:nclass
            optind{curclass} = r.FullInd{curclass}(S{curclass});
        end
        if VERBOSE, fprintf('=> NEW optimum: # Features: %4.0f ==> %s Perf = %g', numel(optind), ActStr, optparam); end
        Opt.S{curclass,cnt} = S; Opt.Param = [Opt.Param optparam];
        cnt=cnt+1;
    else
        if VERBOSE, fprintf('(%s = %g)', ActStr, optparam); end
    end 
    maxK = maxK - numel(rind);
end
         
%% CHECK IF OPTIMIZED FEATURE SPACE PERFORMS BETTER THAN ORIGINAL SPACE
if ~feval(r.evaldir, optparam, r.FullParamMulti)
    optparam = r.FullParam; optind = r.FullInd; optfound = 0; optmodel = r.FullModel;
elseif isnan(param) && ~optfound 
    optind = r.FullInd; optmodel = r.FullModel; optparam = r.FullParamMulti;
    fprintf('\n');warning('Greedy backward search did not return any feature mask for given parameter setting. Return original feature space.')
else
    optfound = 1;
    for curclass=1:nclass
        if r.KneePoint,
            kneepoint = knee_pt(Opt.Param);
            if isnan(kneepoint)
                cprintf('err','\nNot enough data points to compute kneepoint. Selecting final feature mask.');
            else
                fprintf('\nSelected kneepoint of optimization curve at wrapper cycle #%g => %s = %g', kneepoint, ActStr, Opt.Param(kneepoint));
            end
            if isnan(kneepoint), kneepoint = numel(Opt.S{curclass}); end
            optind{curclass} = r.FullInd{curclass}(Opt.S{curclass}{kneepoint});
        else
            optind{curclass} = r.FullInd{curclass}(S{curclass});
        end
       [~,optmodel{curclass}] = feval(TRAINFUNC, Y{curclass}(:,optind{curclass}), label{curclass}, 1, Ps{curclass}); 
    end
end

if VERBOSE; fprintf('\nDone. '); end


