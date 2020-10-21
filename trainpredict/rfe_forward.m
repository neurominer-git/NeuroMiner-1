function [optparam, optind, optfound, optmodel] = ...
    rfe_forward(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr)
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

r = rfe_algo_settings(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr);
S = []; Sind = 1:r.kFea; k = r.kFea; 
optfound = 0; optparam = r.optparam;

if VERBOSE
    fprintf('\n-----------------------------')
    fprintf('\nGREEDY FORWARD FEATURE SEARCH')
    fprintf('\n-----------------------------')
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

lstep = 1;
%% Start Wrapper: FORWARD FEATURE SELECTION
if r.lperc, lstep = ceil((numel(Sind)/100)*r.lperc); end

Opt = struct('S',[],'Param',[],'ParamTs',[]);

switch r.WeightSort 
    
    case 1 %% Sorting is done according to CV1 test performance

        while k > r.MinNum

            lc = numel(Sind);
            if ~lc, break; end
            val = zeros(lc,1); 

            if VERBOSE, fprintf('\n\tFeature pool size: %4.0f out of %4.0f, block size: %4.0f feature(s) ',numel(S), numel(Sind), lstep); end

            while lc > 0

                kS = [S Sind(lc)]; 
                tY = r.Y(:,kS); T = r.T(:,kS);
                [~, model] = feval(TRAINFUNC, tY, label, 1, Ps);
                val(lc) = nk_GetTestPerf(T, r.L, [], model, tY);
                lc = lc - 1;
            end
            
             % Sort features
            [~, ind] = sort(val,r.optfunc);

            % Get lstep features and obtain param for the current feature
            % selection using pre-specified algorithm
            krem = 1:lstep;rind = krem;
            kS = [S Sind(ind(krem))]; 
            tY = r.Y(:,kS); T = r.T(:,kS);
            [~, model] = feval(TRAINFUNC, tY, label, 1, Ps);
            param = nk_GetTestPerf(T, r.L, [], model, tY);
            %param_ts = nk_GetTestPerf(r.Ynew(:,kS), labelnew, [], model, tY);
            
            % Add feature to feature space only if current performance is better
            % then previous space
            % Finally trained model will have only a performance == optparam
            % if the feature stepping is 1
            if feval(r.evaldir, param, optparam) 
                % If user activated random feature selection then we have
                % to randomly choose x% of the features in the block
                if r.FeatRandPerc
                    rstep = ceil(lstep/100)*r.FeatRandPerc;
                    rind = randperm(lstep,rstep);
                % ... otherwise select entire block
                else
                    rind = krem;
                end
                S = [S Sind(ind(rind))]; 
                optparam = param; optfound = 1; 
                if VERBOSE, fprintf('=> NEW optimum: # Features: %4.0f ==> %s = %g', numel(S), ActStr, optparam); end
                Opt.S{end+1} = S; Opt.Param = [Opt.Param optparam];
                %Opt.ParamTs = [Opt.ParamTs param_ts];
            end

            % Remove selected features from feature pool
            Sind(ind(rind)) = [];

            % Recompute lstep according to current feature pool
            if r.lperc, lstep = ceil((numel(Sind)/100)*r.lperc); end
            k = k - numel(rind);
            
        end
        
    case 2
        
         W = abs(nk_GetPrimalW(r.FullModel));
         W = W/(norm(W,2));
         [~, ind] = sort(W,'descend');
         
         while k > r.MinNum
             
            if VERBOSE, fprintf('\n\tFeature pool size: %4.0f out of %4.0f, block size: %4.0f feature(s) ',numel(S), numel(Sind),lstep); end
            
            % Get lstep features 
            krem = 1:lstep; kSind = Sind(ind(krem));
             
            % Extract lstep feature subspace
            kS = [S kSind]; 
            tY = r.Y(:,kS);  T = r.T(:,kS);
             
            % Train model with ks features
            [~, model] = feval(TRAINFUNC, tY, label, 1, Ps);    
            param = nk_GetTestPerf(T, r.L, [], model, tY);
            
            % Add feature to feature space only if current performance is better
            % then previous space
            if feval(r.evaldir, param, optparam) 
                S = [S kSind]; optparam = param;
                Opt.S{end+1} = S; Opt.Param = [Opt.Param optparam];
                if VERBOSE, fprintf('=> NEW optimum: # Features: %4.0f ==> %s = %g', numel(S), ActStr, optparam); end
            else
                if VERBOSE, fprintf('(%s = %g)', ActStr, optparam); end
            end
            
            % Remove selected features from feature pool
            Sind(ind(krem)) = [];
            
            % Retrain model on remaining feature pool
            [~, kmodel] = feval(TRAINFUNC, r.Y(:,Sind), label, 1, Ps);    
            W = abs(nk_GetPrimalW(kmodel));
            [~, ind] = sort(W,'descend');
         
            % Recompute lstep according to current feature pool
            if r.lperc, lstep = ceil((numel(Sind)/100)*r.lperc); end
            k = k - numel(krem);
             
         end
        
end

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
    try
    optind = r.FullInd(Opt.S{kneepoint});
    catch
        fprintf('problem')
    end
elseif isnan(param) && ~optfound 
    % Optimization returned non-finite performance at given parameter
    % combination. Return original feature space
    optind = r.FullInd;
    fprintf('\n');warning('Greedy forward search did not return any feature mask for given parameter setting. Return original feature space.')
else
    optind = r.FullInd(S);
end
if VERBOSE, fprintf('\nDone. '); end
optfound = 1; [~,optmodel] = feval(TRAINFUNC, Y(:,optind), label, 1, Ps); 
