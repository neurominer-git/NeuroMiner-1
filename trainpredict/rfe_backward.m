function [optparam, optind, optfound, optmodel] = ...
    rfe_backward(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr)
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

r = rfe_algo_settings(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr);

%% Feature block size settings
optfound = 0; optparam = r.FullParam; optind = r.FullInd;

S = 1:r.kFea; k = r.kFea; optmodel = r.FullModel;
if k <= r.MinNum, return; end

Opt = struct('S',[],'Param',[],'ParamTs',[]);

%% Start Wrapper: RECURSIVE FEATURE ELIMINATION
if VERBOSE
    fprintf('\n-----------------------------------')
    fprintf('\nGREEDY BACKWARD FEATURE ELIMINATION')
    fprintf('\n-----------------------------------')
    fprintf('\nOptimization data mode: %s', ActStr)
    fprintf('\nParameter evaluation: %s (%s)', r.evaldir, r.optfunc)
    if r.FeatStepPerc
        fprintf('\nStepping: %g%% of %g features per wrapper cycle.',r.lperc, r.kFea)
    else
        fprintf('\nStepping: Top feature in wrapper cycle.')
    end
    if r.FeatRandPerc
        fprintf('\nRandom feature selection: %g%% of top-ranked features in block',r.FeatRandPerc)
    end
end

lstep = ceil((numel(S)/100)*r.lperc);

switch r.WeightSort 
    
    case 1 %% Sorting is done according to CV1 test performance

        while k > r.MinNum 

            l = k;
            lvec = 1:l; lc = numel(lvec);
            if ~lc, break; end
            val = zeros(lc,1); 
            if VERBOSE, fprintf('\n\tFeature pool size: %4.0f, block size: %4.0f feature(s) ',numel(S),lstep); end

            while lc > 0

                lind = lvec; lind(lc)=[]; kS = S(lind);
                tY = r.Y(:,kS); T = r.T(:,kS);
                [~, model] = feval(TRAINFUNC, tY, label, 1, Ps);          
                val(lc) = nk_GetTestPerf(T, r.L, [], model, tY);

                lc = lc - 1;
            end

            % Eliminate feature, that maximally degrades CV1-Test performance
            [param, ind] = sort(val,r.optfunc);
            krem = 1:lstep;
            param = mean(param(krem));
            % If user activated random feature selection then we have
            % to randomly choose x% of the features in the block
            if r.FeatRandPerc
                rstep = ceil(lstep/100)*r.FeatRandPerc;
                rind = randperm(lstep,rstep);
            % ... otherwise select entire block
            else
                rind = krem;
            end
            S(lvec(ind(rind))) = [];
            lstep = ceil((numel(S)/100)*r.lperc);

            if feval(r.evaldir, param, optparam)
                optparam = param; optind = r.FullInd(S); optfound = 1; 
                if VERBOSE, fprintf('=> NEW optimum: # Features: %4.0f ==> %s Perf = %g', numel(optind), ActStr, optparam); end
                Opt.S{end+1} = S; Opt.Param = [Opt.Param optparam];
            else
                    if VERBOSE, fprintf('(%s = %g)', ActStr, optparam); end
                %end
            end 
            k = k - numel(rind);
        end
        
    case 2
        
         W = abs(nk_GetPrimalW(FullModel));
         W = W/(norm(W,2));
         [~, ind] = sort(W,'ascend');
         
         while k > r.MinNum
             
            if VERBOSE, fprintf('\n\tFeature pool size: %g out of %g, block size: %g feature(s) ',numel(S), numel(FullInd), lstep); end
            
            % Get lstep features 
            krem = 1:lstep; 
            
            % Extract lstep feature subspace
            S(ind(krem)) = [];
            tY = r.Y(:,S);  T = r.T(:,S);
             
            % Train and test model with S - krem features
            [~, model] = feval(TRAINFUNC, tY, label, 1, Ps);    
            param = nk_GetTestPerf(T, r.L, [], model, tY);
            
            % Add feature to feature space only if current performance is better
            % then previous space
            if feval(evaldir, param, optparam) 
                optparam = param; optfound = 1;
                if VERBOSE, fprintf('=> NEW optimum: # Features: %4.0f ==> %s = %g', numel(S), ActStr, optparam); end
                Opt.S{end+1} = S; Opt.Param = [Opt.Param optparam];
            end
            
            % Resort features according to current weight vector
            W = abs(nk_GetPrimalW(model));
            [~, ind] = sort(W,'ascend');
         
            % Recompute lstep according to current feature pool
            if r.lperc, lstep = ceil((numel(S)/100)*r.lperc); end
            k = k - numel(krem);
             
         end
    
end

%% CHECK IF OPTIMIZED FEATURE SPACE PERFORMS BETTER THAN ORIGINAL SPACE
if ~feval(r.evaldir, optparam, r.FullParam)
    optparam = r.FullParam; optind = r.FullInd; optfound = 0; optmodel = r.FullModel;
elseif isnan(param) && ~optfound 
    optind = r.FullInd; optmodel = r.FullModel; optparam = r.FullParam;
    fprintf('\n');warning('Greedy backward search did not return any feature mask for given parameter setting. Return original feature space.')
else
    optfound = 1;
    if r.KneePoint,
        kneepoint = knee_pt(Opt.Param,[],true);
        if isnan(kneepoint)
            cprintf('err','\nNot enough data points to compute kneepoint. Selecting final feature mask.');
        else
            fprintf('\nSelected kneepoint of optimization curve at wrapper cycle #%g => %s = %g', kneepoint, ActStr, Opt.Param(kneepoint));
        end
        if isnan(kneepoint), kneepoint = numel(Opt.S); end
        try
            optind = r.FullInd(Opt.S{kneepoint});
        catch
            cprintf('err','\nNo optimum found. Selecting original feature mask.');
            optind = r.FullInd;
        end
    else
        optind = r.FullInd(S);
    end
    [~,optmodel] = feval(TRAINFUNC, Y(:,optind), label, 1, Ps); 
end

if VERBOSE; fprintf('\nDone. '); end


