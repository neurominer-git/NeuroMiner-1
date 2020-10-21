function [opt_hE, opt_E, opt_F, opt_D] = nk_CVMax(E, L, EnsStrat)

global BATCH MODEFL VERBOSE

ind0 = L~=0; E = E(ind0,:); L=L(ind0);
[m,k] = size(E);

% Compute initial ensemble performance
if EnsStrat.Metric == 2, T = sign(E); else T = E; end
opt_hE = nk_EnsPerf(E, L); orig_hE = opt_hE;
opt_D = nk_Entropy(T,[-1 1], m, k); orig_D = opt_D; 
orig_F = 1:k;
if strcmp(MODEFL,'classification')
    varflag = 0;
else
    varflag = 1;
end

switch EnsStrat.type 
    
    case 2
        if ~varflag
            strhdr = 'Binary RCE => max(perf|entropy)';
        else
            strhdr = ['Regression RCE => ' EnsStrat.CompFunc ' (perf|variance)'];
        end
        MaxParam = opt_hE; MaxAbParam = opt_D; I = orig_F; opt_I = I; iter=0;

        % Perform backward elimination that maximizes ambiguity among classifiers
        while k > EnsStrat.MinNum
            l   = k;
            lhE = zeros(l,1);
            leE = zeros(m,l);
            lD  = zeros(l,1);

            while l > 0

                % Compute current ensemble performance and store it in reverse
                % order as "min" and "max" will return the indices of the first
                % element in the vector. However, the last index in our case points
                % to the feature that is least important as defined by
                % nk_CreateSubSets
                kI = I; kI(l) = []; lE = E(:,kI);  lT = T(:,kI); pos = k-l+1;
                % Compute current ensemble performance
                [lhE(pos),leE(:,pos)] = nk_EnsPerf(lE, L);
                % Compute current ensemble ambiguity:
                % in case of classification:    entropy
                % in case of regression:        deviation from mean
                % ensemble expectation
                if ~varflag
                    lD(pos) = nk_Entropy(lT,[-1 1], m, k-1);
                else
                    lD(pos) = nk_RegAmbig(lE, L);
                end
                l=l-1;
            end

            % Find the classifier to be eliminated
            % Check whether performance array consists of identical values
            if length(unique(lhE)) == 1
                % In this case choose the classifier, whose removal increases the
                % ensemble diversity to the maximum degree.
                [abparam, ind] = feval(EnsStrat.CompFunc, lD);
                param = lhE(ind);
            else
                [param,ind] = feval(EnsStrat.CompFunc, lhE);
                abparam = lD(ind);
            end
            
            % Check if elimination of current base learner will
            % lead to complete elimination of this dichotomization
            % in ensemble
            I(k-ind+1) = [];

            if feval(EnsStrat.OptInlineFunc1, param, MaxParam, abparam, MaxAbParam)
                iter        = iter+1;
                opt_I       = I;
                MaxParam    = param;
                MaxAbParam  = abparam;
            end
            k=k-1;
        end
        
    case 6
        
        if ~varflag
            strhdr = 'Binary FCC => max(perf|entropy)';
        else
            strhdr = ['Regression FCC => ' EnsStrat.CompFunc ' (perf|variance)'];
        end
        
        % Create empty set & initialize params for forward search
        origI = 1:k; I=[]; iter=0; opt_I = I;
        switch EnsStrat.CompFunc
            case 'max'
                MaxParam = 0; MaxAbParam = 0; 
            case 'min'
                MaxParam = realmax; MaxAbParam = 0; 
        end
        
        while k > 0
            
            l   = k;
            lhE = zeros(l,1);
            leE = zeros(m,l);
            lD  = zeros(l,1);
            
            while l > 0
                
                pos = k-l+1;
                kI = [I origI(pos)]; lE = E(:,kI);  lT = T(:,kI);
                [lhE(pos),leE(:,pos)] = nk_EnsPerf(lE, L);

                % Compute current ensemble ambiguity (entropy)
                % variance is used as ambiguity criterion in case of
                % regression models
                
                if ~varflag
                    lD(pos) = nk_Entropy(lT,[-1 1], m);
                else
                    lD(pos) = nk_RegAmbig(E(:,I), E(:,origI(pos)));
                end
                l=l-1;
            end
            
            if length(unique(lhE)) == 1
                % In this case choose the classifier, whose removal increases the
                % ensemble diversity to the maximum degree.
                [abparam, ind] = max(lD);
                param = lhE(ind);
            else
                indE = 1:numel(lhE);
                switch EnsStrat.CompFunc
                    case 'min'
                        rlhE = rank_data(lhE,'ascend')'; 
                    case 'max'
                        rlhE = rank_data(lhE,'descend')';
                end
                rlD = rank_data(lD,'descend')';
                rlhElD = mean([rlhE rlD],2); 
                [~,ind] = min(rlhElD); 
                param = lhE(ind);
                %[param,ind] = feval(EnsStrat.CompFunc, lhE);
                abparam = lD(ind);
            end
            
            if feval(EnsStrat.OptInlineFunc1, param, MaxParam, abparam, MaxAbParam)
                iter        = iter+1;
                I           = [I origI(ind)]; % add to the existing feature set
                MaxParam    = param;
                MaxAbParam  = abparam;
                opt_I       = I;
            end
            origI(ind) = [];
            k=k-1;
            
        end 
end
        
if feval(EnsStrat.OptInlineFunc2, opt_hE, MaxParam, opt_D, MaxAbParam)
    rej_str = ' [ solution rejected ]';
    opt_I   = orig_F;
    opt_F   = orig_F;
    opt_E   = E;
else
    rej_str = '';
    opt_hE  = MaxParam;
    opt_D   = MaxAbParam;
    opt_E   = E(:,opt_I);
    opt_F   = orig_F(:,opt_I);
end

if VERBOSE
   fprintf(['\n%s: %g iters\t' ...
    'Ent(orig->final): %1.2f->%1.2f, ' ...
    'Perf (orig->final): %1.2f->%1.2f, ' ...
    '# Learners (orig->final): %g->%g%s'], ...
    strhdr, iter, orig_D, opt_D, orig_hE, opt_hE, numel(orig_F), numel(opt_I), rej_str)
end
