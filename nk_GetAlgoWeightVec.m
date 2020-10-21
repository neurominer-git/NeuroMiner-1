function [xV, AnalP] = nk_GetAlgoWeightVec(SVM, Y, L, MD, decompfl, errorflag)
global STACKING

AnalP = []; xV=[];
switch SVM.prog
    case {'LIBSVM','LIBLIN'}
        switch SVM.kernel.kernstr
            case{' -t 0',' -t 4',' -t 5','lin', 'linear', 'none'} 
                %%%%%%%%%%%%% USE WEIGHTS OF MODEL %%%%%%%%%%%%
                xV = nk_GetPrimalW(MD); % Get weight vector over feature space
                if ~decompfl(1)
                    % Remove cases which are completely NaN
                    [Yn, Ln] = nk_ManageNanCases(Y, L); 
                    AnalP = compute_analytical_pvals(Ln,Yn,xV'); 
                end
            otherwise % non-linear
                %%%%%%%%%% COMPUTE MIN. DIFF. VECTORS %%%%%%%%% 
                xV = nk_VisSV(MD, Y, L);
        end
    case 'MEXELM'
        xV = MD.outW;
    case  'GLMFIT'
        xV = MD.beta(2:end);
    case 'matLRN'
        % Check whether addBias == 1 and whether algo was kernalized
        offs=1; if isfield(MD,'addBias') && MD.addBias, offs = 2; end
        if isfield(MD,'kernel') && MD.kernel, 
            error('Unfortunately, the pre-image of the kernel weight vector cannot be computed with the current version NM.'); 
        end
        if isfield(MD,'w')
            w = MD.w;
        elseif isfield(MD,'weights')
            w = MD.weights;
        else
            error(sprintf('\nNo weight vector found in matLearn model. Unfortunately, the visualization of this model is not supported by NM.'));
        end
        xV = w(offs:end);  
    case 'GLMNET'
        xV = MD.beta(:,end);
    case 'RNDFOR'
        xV = MD.importance;
    case 'DECTRE'
        xV = MD.predictorImportance;
    case 'SEQOPT'
        nF = numel(STACKING.sel_anal);
        xV = zeros(nF,1);
        xV(MD.AnalSeq) = (MD.examsfreq(MD.examsfreq>0)/100)';
    case 'WBLCOX'
        xV = MD.beta;
    case 'IMRELF'
        xV = MD.Weight;
    otherwise
        if errorflag
            error(['Vizualisation for ' SVM.prog ' not supported yet. Please consult the NM Manual for more information']);
        end
end