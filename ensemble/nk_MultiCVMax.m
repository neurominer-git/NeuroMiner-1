% =========================================================================
% FORMAT [opt_hE, opt_E, opt_F, opt_Fcat, opt_D, opt_Pred] = ...
%                                           nk_MultiCVMax(E, L, EnsStrat, C)
% =========================================================================
% 
% Inputs
% ------
% E        :       [m x n] Multi-group Ensemble Hypotheses / Decision Values / 
%                  Probabilities coded in a dichotomization matrix
% L        :       [m x 1] Multi-group Label Vector
% EnsStrat :       Ensemble Strategy settings
% C        :       [1 x n] Dichotomization Vector with 1 => c classes
%
% Outputs 
% -------
% opt_hE   :       Optimized multi-group classification performance
% opt_E    :       Optimized multi-group ensemble
% opt_F    :       Optimized dichotomization vector (cell array with c
%                  cells)
% opt_Fcat :       Optimized dichotomization vector ( optimized C )
% opt_D    :       Optimized Entropy
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011
function [opt_hE, opt_E, opt_F, opt_Fcat, opt_D, opt_Pred] = nk_MultiCVMax(E, L, EnsStrat, C, G)
global BATCH VERBOSE

%% Preparations
ind0    = L~=0; 
E       = E(ind0,:); 
L       = L(ind0);
nclass  = max(C); 

if EnsStrat.Metric == 2, T = sign(E); else T = E; end  

%% Compute initial ensemble performance
[opt_hE, opt_Pred] = nk_MultiEnsPerf(E, T, L, C, G); 
orig_hE = opt_hE;

%% Compute initial ensemble entropy
opt_Dc = zeros(nclass,1);
for curclass=1:nclass
    opt_Dc(curclass) = nk_Entropy(T(:,C==curclass),[1 -1]); 
end
[m, d]  = size(E);
opt_D   = mean(opt_Dc);
orig_D  = opt_D;
opt_F   = cell(nclass,1); 
orig_F  = 1:d; opt_I = orig_F;
k       = zeros(nclass,1); 

for curclass=1:nclass
    k(curclass) = size(E(:,C==curclass),2);
end

ksum = d; 

switch EnsStrat.ConstructMode 

    case 1
        
        %% RECURSIVE FEATURE ELIMINATION
        MaxParam = opt_hE; MaxAbParam = opt_D; iter = 0;
        strhdr = 'Multi-group RCE => max(perf|entropy)';
        I = 1:ksum; MinNum = EnsStrat.MinNum*nclass;
        
        % Perform backward elimination that maximizes ambiguity among classifiers
        while ksum > MinNum

            for curclass = 1:nclass

                l   = k(curclass);
                if ~l, continue; end
                lhE = zeros(l,1);
                lD  = zeros(l,1);
                leE = zeros(m,l);
                indCur = find(C(I)==curclass);

                while l > 0

                    % Compute current ensemble performance and store it in reverse
                    % order as "min" and "max" will return the indices of the first
                    % element in the vector. However, the last index in our case points
                    % to the feature that is least important as defined by
                    % nk_CreateSubSets
                    kX = indCur; kX(l) = []; kI = I; kI(indCur(l)) = []; 
                    lE = E(:,kI); lT = T(:,kI); lC = C(:,kI);
                    pos = k(curclass)-l+1;
                    % Compute current multi-group ensemble performance
                    [lhE(pos), leE(:,pos)] = nk_MultiEnsPerf(lE, lT, L, lC, G);
                    % Compute current dichotomizer ensemble ambiguity (entropy)
                    lD(pos) = nk_Entropy(T(:,kI),[1 -1], m);
                    l=l-1;
                end

                % Find the classifier to be eliminated
                % Check whether performance array consists of identical values
                if length(unique(lhE)) == 1
                    % In this case choose the classifier, whose removal increases the
                    % ensemble diversity to the maximum degree.
                    [abparam, ind]  = max(lD);
                    param           = lhE(ind);
                else
                    [param,ind]     = max(lhE);
                    abparam         = lD(ind);
                end
                topt_Dc             = opt_Dc;
                topt_Dc(curclass)   = abparam;
                topt_D              = mean(topt_Dc);
                
                % Update index vector of base learners
                I(indCur(k(curclass)-ind+1)) = [];

                if feval(EnsStrat.OptInlineFunc1, param, MaxParam, topt_D, MaxAbParam) 
                    iter            = iter+1;
                    opt_Dc          = topt_Dc;
                    MaxParam        = param;
                    MaxAbParam      = topt_D;
                    opt_I           = I;
                end
                ksum = ksum - 1;
                k(curclass) = k(curclass)-1;
            end
        end
            
    case 2
        %% FORWARD CLASSIFIER CONSTRUCTION
        MaxParam = 0; MaxAbParam = 0; MinNum = EnsStrat.MinNum;
        
        strhdr = 'Multi-group FCC => max(perf|entropy)';
        % Initialize set
        I = zeros(1,nclass); origI = []; iter = nclass;
        for curclass=1:nclass
            ind = find(C==curclass);
            I(curclass) = ind(1);
            origI = [origI ind(2:end)];
            k(curclass) = k(curclass)-1;
        end
        ksum = ksum-nclass;
        
        while ksum > MinNum
            
            for curclass = 1:nclass
            
                l        = k(curclass);
                if ~l, continue; end
                lhE      = zeros(l,1);
                lD       = zeros(l,1);
                indCur   = find(C(origI)==curclass);
                
                while l > 0
                    
                    % Compute current component index
                    kI = [I origI(indCur(l))];
                    lE = E(:,kI); lT = T(:,kI); lC = C(:,kI);
                    
                    % Compute current ensemble performance
                    lhE(l) = nk_MultiEnsPerf(lE, lT, L, lC, G);

                    % Compute current ensemble ambiguity (entropy)
                    lD(l) = nk_Entropy(lT(:,lC==curclass),[-1 1], m);
                    l=l-1;
                end

                if length(unique(lhE)) == 1
                    % In this case choose the classifier, whose addition increases the
                    % ensemble diversity to the maximum degree.
                    [abparam, ind] = max(lD);
                    param = lhE(ind);
                else
                    [param,ind] = max(lhE);
                    abparam = lD(ind);
                end 
                
                topt_Dc             = opt_Dc;
                topt_Dc(curclass)   = abparam;
                topt_D              = mean(topt_Dc);
                
                if feval(EnsStrat.OptInlineFunc1, param, MaxParam, topt_D, MaxAbParam) 
                    iter        = iter+1;
                    I           = [I origI(indCur(ind))]; % add to the existing feature set
                    MaxParam    = param;
                    MaxAbParam  = topt_D;
                    opt_I       = I;
                end
                
                origI(indCur(ind)) = [];
                k(curclass) = k(curclass)-1;
                ksum = ksum - 1 ;
            end
            
        end
end

% Check if optimized ensemble achieved better results than original
% ensemble: This can be a) a better classification performance, or b) an
% equal classification performance but a better diversity measure

if feval(EnsStrat.OptInlineFunc2, orig_hE, MaxParam, orig_D, MaxAbParam)
    opt_I = orig_F; opt_hE = orig_hE; opt_D = orig_D;
else
    opt_D = MaxAbParam; 
    [opt_hE, opt_Pred] = nk_MultiEnsPerf(E(:,opt_I), T(:,opt_I), L, C(:,opt_I), G); 
end

% Generate outputs of function
opt_E = E(:,opt_I); log_I = false(1,d); log_I(opt_I) = true; opt_Fcat = opt_I;
for curclass=1:nclass
    indCur = C == curclass; opt_F{curclass} = log_I(indCur);
end

% Print some info
if VERBOSE
    fprintf(['\n%s: %g iters\t' ...
        'Ent(orig->final): %1.2f->%1.2f, ' ...
        'Perf (orig->final): %1.2f->%1.2f, ' ...
        '# Learners (orig->final): %g->%g'], ...
        strhdr, iter, orig_D, opt_D, orig_hE, opt_hE, numel(orig_F), numel(opt_I))
end

