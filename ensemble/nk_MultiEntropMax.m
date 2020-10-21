function [opt_hE, opt_E, opt_F, opt_Fcat, opt_D] = nk_MultiEntropMax(E, L, EnsStrat, C, G)
global BATCH

ind0    = L~=0; 
E       = E(ind0,:); 
L       = L(ind0);
nclass  = max(C); 

if EnsStrat.Metric == 2, T = sign(E); else T = E; end

% Compute initial ensemble performance
opt_hE  = nk_MultiEnsPerf(E, T, L, C, G); 
orig_hE = opt_hE;

% Compute initial ensemble ambiguity
opt_Dc = zeros(nclass,1);
for curclass=1:nclass
    opt_Dc(curclass) = nk_Lobag(T(:,C==curclass)); 
end
[m, d]  = size(E);
opt_D   = mean(opt_Dc);
orig_D  = opt_D;
opt_F   = cell(nclass,1); 
orig_F  = 1:d; opt_I = orig_F;
iter    = 0;
k       = zeros(nclass,1); 

switch EnsStrat.type 
    case 1
        MaxParam = opt_D; 
        MinNum = EnsStrat.MinNum*nclass;
    case 5
        MaxParam = 0; 
end

for curclass=1:nclass
    k(curclass) = size(E(:,C==curclass),2);
end

ksum = d;
switch EnsStrat.type

    case 1

        strhdr = 'Multi-group RCE => max(entropy)';
        I = 1:ksum;
        
        % Perform backward elimination that maximizes ambiguity among classifiers
        while ksum > MinNum

            for curclass = 1:nclass

                l   = k(curclass);
                lD  = zeros(l,1);
                indCur = find(C(I)==curclass);

                while l > 0

                    % Compute current ensemble performance and store it in reverse
                    % order as "min" and "max" will return the indices of the first
                    % element in the vector. However, the last index in our case points
                    % to the feature that is least important as defined by
                    % nk_CreateSubSets
                    kX = indCur; kX(l) = []; pos = k(curclass)-l+1;
                    
                    % Compute current dichotomizer ensemble ambiguity (entropy)
                    lD(pos) = nk_Entropy(T(:,kX),[1 -1], m);
                    l=l-1;
                end

                % Find the classifier to be eliminated
                [param, ind]  = max(lD);
                topt_Dc             = opt_Dc;
                topt_Dc(curclass)   = param;
                topt_D              = mean(topt_Dc);
                
                % Update index vector of base learners
                I(indCur(k(curclass)-ind+1)) = [];

                if topt_D >= MaxParam
                    iter            = iter+1;
                    opt_Dc          = topt_Dc;
                    MaxParam        = topt_D;
                    opt_I           = I;
                end
                ksum = ksum - 1;
                k(curclass) = k(curclass)-1;
            end
        end
            
    case 5

        strhdr = 'Multi-group FCC => max(entropy)';
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
                    
                    kI = [I origI(indCur(l))]; lT = T(:,kI); lC = C(:,kI);
                    % Compute current ensemble ambiguity (entropy)
                    lD(l) = nk_Entropy(lT(:,lC==curclass),[-1 1], m);
                    l=l-1;
                end

                [param, ind] = max(lD);
                topt_Dc             = opt_Dc;
                topt_Dc(curclass)   = param;
                topt_D              = mean(topt_Dc);
                
                if topt_D > MaxParam
                    iter        = iter+1;
                    I           = [I origI(indCur(ind))]; % add to the existing feature set
                    MaxParam    = topt_D;
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
if opt_D >= MaxParam
    opt_I = orig_F; 
else
    opt_D = MaxParam;
    opt_hE  = nk_MultiEnsPerf(E(:,opt_I), T(:,opt_I), L, C(:,opt_I), G); 
end

% Generate outputs
opt_E = E(:,opt_I); log_I = false(1,d); log_I(opt_I) = true; opt_Fcat = opt_I;
for curclass=1:nclass
    indCur = C == curclass; opt_F{curclass} = log_I(indCur);
end

% Print some info
if ~BATCH
fprintf(['\n%s: %g iters\t' ...
    'Ent(orig->final): %1.2f->%1.2f, ' ...
    'Perf (orig->final): %1.2f->%1.2f, ' ...
    '# Learners (orig->final): %g->%g'], ...
    strhdr, iter, orig_D, opt_D, orig_hE, opt_hE, numel(orig_F), numel(opt_I))
end

return