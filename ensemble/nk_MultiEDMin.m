% =========================================================================
% FORMAT [opt_hE, opt_E, opt_F, opt_Fcat, opt_D, opt_Pred] = ...
%                                           nk_MultiEDMin(E, L, EnsStrat, C)
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
% opt_D    :       Optimized EDmin ( ED = bias - unbiased var + biased var )
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [opt_hE, opt_E, opt_F, opt_Fcat, opt_D, opt_Pred] = nk_MultiEDMin(E, L, EnsStrat, C, G)
global BATCH

ind0    = L~=0; 
E       = E(ind0,:); 
L       = L(ind0);
nclass  = max(C); 
uC      = unique(C);
if EnsStrat.Metric == 2, T = sign(E); else T = E; end

% Compute initial ensemble performance
[opt_hE, opt_Pred]  = nk_MultiEnsPerf(E, T, L, C, G); 
orig_hE = opt_hE;

% Compute initial ensemble ambiguity
opt_D = nk_LobagMulti(E, T, L, C, uC, nclass, G);

[m, d]  = size(E);
orig_D  = opt_D;
opt_F   = cell(nclass,1); 
orig_F  = 1:d; opt_I = orig_F;
iter    = 0;
k       = zeros(nclass,1); 

switch EnsStrat.ConstructMode 
    case 1
        MaxParam = opt_D; 
        MinNum = EnsStrat.MinNum*nclass;
    case 2
        MaxParam = Inf; 
end

for curclass=1:nclass
    k(curclass) = size(E(:,C == curclass),2);
end

ksum = d;
switch EnsStrat.ConstructMode 

    case 1

        strhdr = 'Multi-group RCE => min(bias+variance)';
        I = 1:ksum;
        
        % Perform backward elimination that maximizes ambiguity among classifiers
        while ksum > MinNum

            for curclass = 1:nclass

                l   = k(curclass);
                lD  = zeros(l,1);
                indCur = find(C(I)==curclass);
                %indNonCur = find(C(I) ~= curclass);
                
                while l > 0

                    % Compute current ensemble performance and store it in reverse
                    % order as "min" and "max" will return the indices of the first
                    % element in the vector. However, the last index in our case points
                    % to the feature that is least important as defined by
                    % nk_CreateSubSets
                    kI = I; kI(indCur(l)) = []; 
                    %kX = indCur; kX(l) = []; kX = [kX indNonCur]; 
                    
                    lD( k(curclass)-l+1 ) = nk_LobagMulti(E(:,kI),T(:,kI), L, C(kI), uC, nclass, G);
                    l=l-1;
                end

                % Find the classifier to be eliminated
                [param, ind]        = min(lD);
               
                % Update index vector of base learners
                I(indCur(k(curclass)-ind+1)) = [];

                if param <= MaxParam
                    iter            = iter+1;
                    MaxParam        = param;
                    opt_I           = I;
                end
                ksum = ksum - 1;
                k(curclass) = k(curclass)-1;
            end
        end
            
    case 2

        strhdr = 'Multi-group FCC => min(bias+variance)';
        % Initialize set
        I = zeros(1,nclass); origI = [];
        for curclass=1:nclass
            ind = find(C==curclass);
            I(curclass) = ind(1);
            origI = [origI ind(2:end)];
            k(curclass) = k(curclass)-1;
        end
        ksum = ksum-nclass;
        
        while ksum > 0
            
            for curclass = 1:nclass
            
                l        = k(curclass);
                lD       = zeros(l,1);
                indCur   = find(C(origI)==curclass);
                
                while l > 0
                    
                    kI = [I origI(indCur(l))];
                    lD(l) = nk_LobagMulti(E(:,kI), T(:,kI), L, C(kI), uC, nclass, G);
                    l=l-1;
                end

                [param, ind] = min(lD);
                
                if param < MaxParam
                    iter        = iter+1;
                    I           = [I origI(indCur(ind))]; % add to the existing feature set
                    MaxParam    = param;
                    opt_I       = I;
                end
                origI(indCur(ind))  = [];
                k(curclass) = k(curclass)-1;
                ksum        = ksum - 1;
            end
        end
end

% Check if optimized ensemble achieved better results than original
% ensemble
if orig_D < MaxParam
    opt_I = orig_F; 
else
    opt_D = MaxParam;
    [opt_hE, opt_Pred] = nk_MultiEnsPerf(E(:,opt_I), T(:,opt_I), L, C(:,opt_I), G); 
end

% Generate outputs
opt_E = E(:,opt_I); log_I = false(1,d); log_I(opt_I) = true; opt_Fcat = opt_I;
for curclass=1:nclass
    indCur = C == curclass; opt_F{curclass} = log_I(indCur);
end

% Print some info
if ~BATCH
fprintf(['\n%s: %g iters\t' ...
    'Div(orig->final): %1.2f->%1.2f, ' ...
    'Perf (orig->final): %1.2f->%1.2f, ' ...
    '# Learners (orig->final): %g->%g'], ...
    strhdr, iter, orig_D, opt_D, orig_hE, opt_hE, numel(orig_F), numel(opt_I))
end

return
