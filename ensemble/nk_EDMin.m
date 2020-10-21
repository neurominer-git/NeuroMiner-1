function [opt_hE, opt_E, opt_F, opt_D] = nk_EDMin(E, L, EnsStrat)

global BATCH

ind0 = L~=0; E = E(ind0,:); L=L(ind0);
[m,k] = size(E);

% Compute initial ensemble performance
if EnsStrat.Metric == 2, T = sign(E); else T = E; end
opt_hE = nk_EnsPerf(E, L); orig_hE = opt_hE;
opt_D = nk_Lobag(T, L); orig_D = opt_D; 
orig_F = 1:k;

switch EnsStrat.type 
    
    case 4
        
        strhdr = 'Binary RCE => min(lobag)';
        EvalParam = opt_D; I = orig_F; opt_I = I; iter=0;

        % Perform backward elimination that maximizes ambiguity among classifiers
        while k > EnsStrat.MinNum
            
            l = k; lD = zeros(l,1);

            while l > 0

                % Compute current ensemble performance and store it in reverse
                % order as "min" and "max" will return the indices of the first
                % element in the vector. However, the last index in our case points
                % to the feature that is least important as defined by
                % nk_CreateSubSets
                kI = I; kI(l) = []; lT = T(:,kI); pos = k-l+1;
               
                % Compute current ensemble ambiguity (diversity)
                lD(pos) = nk_Lobag(lT, L);
                
                l=l-1;
            end

            % Find the classifier to be eliminated
            % Choose the classifier, whose removal increases the
            % ensemble diversity to the maximum degree.
            [param, ind] = min(lD);
                
            I(k-ind+1) = [];

            if param <= EvalParam
                iter        = iter+1;
                opt_I       = I;
                EvalParam   = param;
            end
            k=k-1;
        end
        
    case 8
        
        strhdr = 'Binary FCC => min(lobag)';
        
        % Create empty set
        origI = 1:k; I=[]; EvalParam = Inf; iter=0;
        
        while k > 0
            
            l = k; lD = zeros(l,1);
            
            while l > 0
                
                pos = k-l+1;
                kI = [I origI(pos)]; lT = T(:,kI);
                
                % Compute current ensemble ambiguity (entropy)
                lD(pos) = nk_Lobag(lT, L);
                l=l-1;
            end
            
            [param,ind] = min(lD);
            
            if param <= EvalParam
                iter        = iter+1;
                I           = [I origI(ind)]; % add to the existing feature set
                EvalParam    = param;
                opt_I       = I;
            end
            origI(ind) = [];
            k=k-1;
        end 
end
        
if opt_D < EvalParam 
    opt_I   = orig_F; 
    opt_F   = orig_F;
    opt_E   = E;
else
    opt_E   = E(:,opt_I);
    opt_hE  = nk_EnsPerf(opt_E, L);
    opt_D   = EvalParam;
    opt_F   = orig_F(opt_I);
end

if ~BATCH
   fprintf(['\n%s: %g iters\t' ...
    'Div(orig->final): %1.2f->%1.2f, ' ...
    'Perf (orig->final): %1.2f->%1.2f, ' ...
    '# Learners (orig->final): %g->%g'], ...
    strhdr, iter, orig_D, opt_D, orig_hE, opt_hE, numel(orig_F), numel(opt_I))
end


return