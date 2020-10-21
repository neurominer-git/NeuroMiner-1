function nk_AutomatePermLeaveSiteOutSamples(tNM, act, M, analdim)
% =========================================================================
% AUTOMATION EXAMPLE 2:
% tNM :     template NM structure
% act :     action mode:
%           'create_M' => create permutation index array (see below)
%           'create_NM' => update template NM according to permutations and
%           save updated NM
%           'run_NM' => run permuted NM analyses
%           'get_results_NM' => retrieve results from permuted NM analyses
% M :       permutation index array:
%           dim 1: cases: 1 = move to OOCV,0 = retain in training/CV data, 
%           dim 2: permutations of dim 1, 
%           dim 3: different sizes of subsamples to be moved to OOCV)
% analdim : Index to template analysis
% -------------------------------------------------------------------------
% Use M (permutation index array) to compare N random leave-group-out 
% partitions against originally observed leave-group-out CV. The original 
% leave-group-out setting should have been defined in the template NM
% (=>tNM) structure (NM.TrainParam.RAND) and should - of course - be present 
% in the given analysis containing the performance measures derived from the
% observed leave-group-out CV (=>analdim).
% If M has not been defined yet invoke the script with act = 'create_M'.
% Then, load the Permutation index array from the resulting file and rerun
% the script with 'create_NM', 'run_NM' or 'get_results_NM'.
% After defining M, the automation script (create_NM) recreates the CV 
% structure. The NM structure updated on each permutation vector of M is 
% saved to disk. Then, the script can be invoked with the 'run_NM' action 
% mode to compute the analyses. Finally, 'get_results_NM' can be used to 
% integrate the analysis results across all computed NM models.
% =========================================================================
% (c) N. Koutsouleris V0.1, 05/2017
global NM BATCH

groups = []; appendfl=false; oldcv=[];

if strcmp(act,'create_M')
    if exist('M','var') && ~isempty(M)
        warning('Recreating permutation index array')
    end
    L = size(tNM.label,1);
    O = nk_input('Define original site membership vector','e',[],L);
    uO = nk_CountUniqueOccurences(O); fprintf('\nFound %g sites in vector',numel(uO));
    lN = nk_input('Define number of random permutations of the site vector',0,'e',1000);
    M=[];
    for i=1:lN
        I = randperm(L); 
        M = [M O(I)]; 
    end 
    save('Perm_Index_Array.mat','M')
    return
end

if ~exist('M','var') || isempty(M)
    error('Permutation index array needed.')
end

if isfield(tNM,'analysis')
    if ~exist('analdim','var') || isempty(analdim)
        n_anal = length(tNM.analysis);
        prmpt = sprintf('Found %g analyses in template NM. Indicate on which analysis to work on!',n_anal);
        analdim = nk_input(prmpt,0,'e',[],[1 n_anal]);
    end
else
    error('No analysis found in NM structure. Setup and initialize an analysis first before invoking automation script!')
end

mM = size(M,2);

if strcmp(act,'get_results_NM')
    TP_CV = zeros(mM,1);
    TN_CV = zeros(mM,1);
    FP_CV = zeros(mM,1);
    FN_CV  = zeros(mM,1);
    Sens_CV = zeros(mM,1);
    Spec_CV = zeros(mM,1);
    BAC_CV  = zeros(mM,1);
    PPV_CV  = zeros(mM,1);
    NPV_CV  = zeros(mM,1);
    PSI_CV  = zeros(mM,1);
    pLR_CV  = zeros(mM,1);
    DOR_CV  = zeros(mM,1);
    NND_CV  = zeros(mM,1);
    uO = nk_CountUniqueOccurences(M(:,1)); nO = numel(uO);
    TP_lCV = zeros(mM,nO);
    FP_lCV = zeros(mM,nO);
    TN_lCV = zeros(mM,nO);
    FN_lCV = zeros(mM,nO);
    Sens_lCV = zeros(mM,nO);
    Spec_lCV = zeros(mM,nO);
    BAC_lCV = zeros(mM,nO);
    PPV_lCV = zeros(mM,nO);
    NPV_lCV = zeros(mM,nO);
    PSI_lCV = zeros(mM,nO);
    pLR_lCV = zeros(mM,nO);
    DOR_lCV = zeros(mM,nO);
    NND_lCV = zeros(mM,nO);
end

for j=1:mM

    fNM = sprintf('NM_Perm%g.mat',j);
    pNM = fullfile(pwd,fNM);
 
    switch act

        case 'create_NM'
            
            if ~exist(pNM,'file')
                fprintf('\nCreating NM permutation file %s',pNM);
                NM = tNM;
                NM.Y{1}=[];
                NM.TrainParam.RAND.CV2LCO.ind = M(:,j);
                NM.TrainParam.RAND.InnerFold = 10;
                NM.cv = nk_MakeCrossFolds(NM.label, NM.TrainParam.RAND, NM.modeflag, groups, NM.groupnames, oldcv, appendfl, 1);
                NM.analysis{analdim}.params.cv = NM.cv;
                NM.analysis{analdim}.params.TrainParam.SAV.matname = sprintf('RESIS-PANSS-NS_NonResp_vs_Resp_Classifier_Perm%g',j);
                save(pNM,'NM');
            else
                fprintf('\nNM permutation file %s exists. Skipping',pNM);
            end

        case 'run_NM'
            
            strout = nk_Preprocess_StrCfg([], []);
            load(pNM);
            CVdimanalysis = fullfile(NM.analysis{analdim}.rootdir, 'LIBSVM', [ NM.analysis{analdim}.params.TrainParam.SAV.matname ...
                                                    '_CVdimanalysis' strout '_ID' NM.id '.mat']); 
            if exist(CVdimanalysis,'file'), 
                fprintf('\n%s exists. Skipping permutation %g.', CVdimanalysis, j)
                continue; 
            end
            BATCH = true;
            NM.Y{1} = tNM.Y{1};
            uO = nk_CountUniqueOccurences(NM.TrainParam.RAND.CV2LCO.ind); fprintf('\nFound %g sites in NM structure',numel(uO));
            NM.TrainParam.RAND.OuterFold = numel(uO);
            NM.TrainParam.RAND.OuterPerm = 1;
            NM.TrainParam.RAND.InnerPerm = 1;
            inp = struct('analind',1, 'lfl',1,'preprocmat',[],'gdmat',[], 'gdanalmat', [], 'varstr', [], 'concatfl', [], 'ovrwrt', 2, 'update', true);
            inp.GridAct = true(NM.TrainParam.RAND.OuterPerm, NM.TrainParam.RAND.OuterFold );
            [~, ~, NM] = nk_MLOptimizerPrep(7, NM, inp, 'automation_script');
            %NM.analysis{analdim} = nk_MLOptimizerPrep(Inf, NM, inp, par);
            save(pNM,'NM');

        case 'get_results_NM'
             fprintf('\nLoading %s',pNM); load(pNM);
             % Compute global performance
             TP_CV(j) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.TP;
             FP_CV(j) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.FP;
             TN_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.TN;
             FN_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.FN;
             Sens_CV(j) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.sens;
             Spec_CV(j) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.spec;
             BAC_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.BAC;
             PPV_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.PPV;
             NPV_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.NPV;
             PSI_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.PSI;
             pLR_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.pLR;
             DOR_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.DOR;
             NND_CV(j)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.NND;
             % Compute simulated leave-site out performances
             [uO, fO] = nk_CountUniqueOccurences(M(:,j)); [uO, s_ind] = sort(uO); fO = fO(s_ind);
             for i=1:numel(uO)
                 ind = M(:,j)==uO(i);
                 Li = NM.label(ind); Li(Li==2)=-1;
                 Pi = NM.analysis{analdim}.GDdims{1}.BinClass{1}.mean_predictions(ind);
                 Ci = ALLPARAM(Li,Pi);
                 TP_lCV(j,i) = Ci.TP;
                 FP_lCV(j,i) = Ci.FP;
                 TN_lCV(j,i)  = Ci.TN;
                 FN_lCV(j,i)  = Ci.FN;
                 Sens_lCV(j,i) = Ci.sens;
                 Spec_lCV(j,i) = Ci.spec;
                 BAC_lCV(j,i)  = Ci.BAC;
                 PPV_lCV(j,i)  = Ci.PPV;
                 NPV_lCV(j,i)  = Ci.NPV;
                 PSI_lCV(j,i)  = Ci.PSI;
                 pLR_lCV(j,i)  = Ci.pLR;
                 DOR_lCV(j,i)  = Ci.DOR;
                 NND_lCV(j,i)  = Ci.NND;
             end
    end
end

if strcmp(act,'get_results_NM')
    save Perm_Results.mat TP_CV FP_CV TN_CV FN_CV Sens_CV Spec_CV BAC_CV PPV_CV NPV_CV PSI_CV pLR_CV DOR_CV NND_CV ...
        TP_lCV FP_lCV TN_lCV FN_lCV Sens_lCV Spec_lCV BAC_lCV PPV_lCV NPV_lCV PSI_lCV pLR_lCV DOR_lCV NND_lCV
    H = {'TP','FP','TN','FN','Sens','Spec','BAC','PPV','NPV' ,'PSI','pLR','DOR','NND'};
    DOR_CV(isinf(DOR_CV)) = NaN; 
    pLR_CV(isinf(pLR_CV)) = NaN;
    NND_CV(isinf(NND_CV)) = NaN;
    DOR_lCV(isinf(DOR_lCV)) = max(DOR_lCV(isfinite(DOR_lCV)));
    pLR_lCV(isinf(pLR_lCV)) = max(pLR_lCV(isfinite(pLR_lCV)));
    NND_lCV(isinf(NND_lCV)) = max(NND_lCV(isfinite(NND_lCV)));
    A = nm_nanmedian(TP_CV);            B = nm_nanstd(TP_CV);
    A = [A; nm_nanmean(FP_CV)];       B = [B; nm_nanstd(FP_CV)];
    A = [A; nm_nanmean(TN_CV)];       B = [B; nm_nanstd(TN_CV)];
    A = [A; nm_nanmean(FN_CV)];       B = [B; nm_nanstd(FN_CV)];
    A = [A; nm_nanmean(Sens_CV)];     B = [B; nm_nanstd(Sens_CV)];
    A = [A; nm_nanmean(Spec_CV)];     B = [B; nm_nanstd(Spec_CV)]; 
    A = [A; nm_nanmean(BAC_CV)];      B = [B; nm_nanstd(BAC_CV)];
    A = [A; nm_nanmean(PPV_CV)];      B = [B; nm_nanstd(PPV_CV)];
    A = [A; nm_nanmean(NPV_CV)];      B = [B; nm_nanstd(NPV_CV)];
    A = [A; nm_nanmean(PSI_CV)];      B = [B; nm_nanstd(PSI_CV)];
    A = [A; nm_nanmean(pLR_CV)];      B = [B; nm_nanstd(pLR_CV)];
    A = [A; nm_nanmean(DOR_CV)];      B = [B; nm_nanstd(DOR_CV)];
    A = [A; nm_nanmean(NND_CV)];      B = [B; nm_nanstd(NND_CV)];
    A = A';                          B = B';
    Mean_CV_tbl = array2table(A,'VariableNames',H);
    SD_CV_tbl = array2table(B,'VariableNames',H);
    writetable(Mean_CV_tbl,'Perm_Results.xls','Sheet','Mean_CV');
    writetable(SD_CV_tbl,'Perm_Results.xls','Sheet','SD_CV');
    for i=1:nO        
        A = nm_nanmedian(TP_lCV(:,i));            B = nm_nanstd(TP_lCV(:,i));
        A = [A; nm_nanmean(FP_lCV(:,i))];       B = [B; nm_nanstd(FP_lCV(:,i))];
        A = [A; nm_nanmean(TN_lCV(:,i))];       B = [B; nm_nanstd(TN_lCV(:,i))];
        A = [A; nm_nanmean(FN_lCV(:,i))];       B = [B; nm_nanstd(FN_lCV(:,i))];
        A = [A; nm_nanmean(Sens_lCV(:,i))];     B = [B; nm_nanstd(Sens_lCV(:,i))];
        A = [A; nm_nanmean(Spec_lCV(:,i))];     B = [B; nm_nanstd(Spec_lCV(:,i))]; 
        A = [A; nm_nanmean(BAC_lCV(:,i))];      B = [B; nm_nanstd(BAC_lCV(:,i))];
        A = [A; nm_nanmean(PPV_lCV(:,i))];      B = [B; nm_nanstd(PPV_lCV(:,i))];
        A = [A; nm_nanmean(NPV_lCV(:,i))];      B = [B; nm_nanstd(NPV_lCV(:,i))];
        A = [A; nm_nanmean(PSI_lCV(:,i))];      B = [B; nm_nanstd(PSI_lCV(:,i))];
        A = [A; nm_nanmean(pLR_lCV(:,i))];      B = [B; nm_nanstd(pLR_lCV(:,i))];
        A = [A; nm_nanmean(DOR_lCV(:,i))];      B = [B; nm_nanstd(DOR_lCV(:,i))];
        A = [A; nm_nanmean(NND_lCV(:,i))];      B = [B; nm_nanstd(NND_lCV(:,i))];
        A = A';                                B = B';
        Mean_CV_tbl = array2table(A,'VariableNames',H);
        SD_CV_tbl = array2table(B,'VariableNames',H);
        writetable(Mean_CV_tbl,'Perm_Results.xls','Sheet',sprintf('Mean_CV_L%g',i) );
        writetable(SD_CV_tbl,'Perm_Results.xls','Sheet',sprintf('SD_CV_L%g',i));
    end
    % Compute significance
    oBAC_CV = tNM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.BAC;
    % Probability for observed BAC (oBAC) being greater or equal permuted
    % BAC => Every time the observed BAC is lower, this increases the
    % probability that the H0 (Leave-Group-Out CV is worse than pooled CV) 
    % is rejected! Do following at the global level (across all sites): Test H1
    Pval_BAC_CV = sum( BAC_CV <= oBAC_CV ) / mM;
    VarNames = cell(nO+1,1); VarNames{1}='All'; RowNames = {'N','P value'};
    % ... and now for each site
    for i=1:nO
        [~,ind] = max(NM.covars(:,1:3),[],2);
        Li = tNM.label(ind==uO(i)); Li(Li==2)=-1;
        Pi = tNM.analysis{analdim}.GDdims{1}.BinClass{1}.mean_predictions(ind==uO(i));
        Ci = ALLPARAM(Li,Pi);
        oBAC_lCV = Ci.BAC;
        Pval_BAC_lCV(i) = sum( BAC_lCV(:,i) <= oBAC_lCV  ) / mM;
        VarNames{i+1} = sprintf('Group_%g',i);
    end
    Ns = [ sum(fO) fO' ];
    Pvals = [Pval_BAC_CV Pval_BAC_lCV];
    tblarr = [Ns; Pvals];
    pval_tbl = array2table(tblarr, 'VariableNames', VarNames, 'RowNames', RowNames);
    writetable(pval_tbl,'Perm_Results.xls','Sheet','Results');
end
clear NM