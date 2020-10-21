function nk_AutomatePermOOCVSamples(tNM, act, M, analdim)
% =========================================================================
% AUTOMATION EXAMPLE 1:
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
% Use M (permutation index array) to move subsample of given template NM (tNM) 
% population to independent test data (OOCV container should have been created 
% in NM already). First, the automation script (create_NM) recreates the CV 
% structure and updates some preprocessing settings (should be adjusted by 
% the user) so that the automation will properly work. The NM structure 
% updated on each permutation vector of M is saved to disk. Then, the
% script can be invoked with the 'run_NM' action mode to compute the
% analyses. Finally, 'get_results_NM' can be used to integrate the analysis
% results across all computed NM models.
% =========================================================================
% (c) N. Koutsouleris V0.1, 12/2016
global NM

O.fldnam = 'OOCV';
O.ind = 1;
O.copycut = 1;
O.mode = 2;
groups = [];
appendfl=false; 
oldcv=[];

if strcmp(act,'create_M')
    if exist('M','var') && ~isempty(M)
        warning('Recreating permutation index array')
    end
    L = size(tNM.label,1);
    T = nk_input('Define vector with N cases to be moved to OOCV',0,'e',[],[1 L]);
    lT = length(T); 
    lN = nk_input('Define number of permutations for each OOCV-N',0,'e',50);
    M = false(L,lN,lT);
    for j=1:lT, for i=1:lN, I = randperm(L,T(j)); M(I,i,j) = true; end; end
    save('Perm_Index_Array.mat','M')
    return
end

if ~exist('M','var') || isempty(M)
    error('Permutation index array needed.')
end

if isfield(NM,'analysis')
    if ~exist('analdim','var') || isempty(analdim)
        n_anal = length(tNM.analysis);
        prmpt = sprintf('Found %g analyses in template NM. Indicate on which analysis to work on!',n_anal);
        analdim = nk_input(prmpt,0,'e',[],[1 n_anal]);
    end
else
    error('No analysis found in NM structure.\nSetup and initialize an analysis first invoking automation script!')
end

nM = size(M,3);
mM = size(M,2);

if strcmp(act,'get_results_NM')
    TP_CV = zeros(mM,nM);
    TN_CV = zeros(mM,nM);
    FP_CV = zeros(mM,nM);
    FN_CV  = zeros(mM,nM);
    Sens_CV = zeros(mM,nM);
    Spec_CV = zeros(mM,nM);
    BAC_CV  = zeros(mM,nM);
    PPV_CV  = zeros(mM,nM);
    NPV_CV  = zeros(mM,nM);
    PSI_CV  = zeros(mM,nM);
    pLR_CV  = zeros(mM,nM);
    DOR_CV  = zeros(mM,nM);
    NND_CV  = zeros(mM,nM);
    
    TP_OOCV = zeros(mM,nM);
    TN_OOCV = zeros(mM,nM);
    FP_OOCV = zeros(mM,nM);
    FN_OOCV  = zeros(mM,nM);
    Sens_OOCV = zeros(mM,nM);
    Spec_OOCV = zeros(mM,nM);
    BAC_OOCV  = zeros(mM,nM);
    PPV_OOCV  = zeros(mM,nM);
    NPV_OOCV  = zeros(mM,nM);
    PSI_OOCV  = zeros(mM,nM);
    pLR_OOCV  = zeros(mM,nM);
    DOR_OOCV  = zeros(mM,nM);
    NND_OOCV  = zeros(mM,nM);
end

for i=1:nM
    
    for j=1:mM
        O.vec = M(:,j,i);
        NVEC = sum(O.vec);
        fNM = sprintf('NM_Perm%g_N%g.mat',j,NVEC);
        pNM = fullfile(pwd,fNM);
        
        switch act
        
            case 'create_NM'
              
                O.desc = sprintf('Perm%g_N%g.mat',j,NVEC);
                O.date = date;
               
                NM = nk_DefineOOCVData_config(tNM, O, 'Permutator');
                [ID, NUM] = nk_CountUniqueOccurences(NM.label);
                minF_CV2 = min(NUM);
                minF_CV1 = min(NUM)-1;
                NM.TrainParam.RAND.OuterFold = minF_CV2;
                NM.TrainParam.RAND.OuterPerm = 1;
                NM.TrainParam.RAND.InnerFold = minF_CV1;
                NM.TrainParam.RAND.InnerPerm = 1;
               
                NM.cv = nk_MakeCrossFolds(NM.label, NM.TrainParam.RAND, NM.modeflag, groups, NM.groupnames, oldcv, appendfl);
                NM.analysis{analdim}.params.cv = NM.cv;
                NM.analysis{analdim}.params.TrainParam.PREPROC{5}.ACTPARAM{2}.REMVARCOMP.dimmode=3;
                NM.analysis{analdim}.params.TrainParam.PREPROC{5}.ACTPARAM{2}.REMVARCOMP.G(M(:,j,i),:)=[];
                NM.analysis{analdim}.params.TrainParam.SAV.matname = sprintf('RESIS-PANSS-NS_NonResp_vs_Resp_Classifier_Perm%g_N%g',j,NVEC);
               
                save(pNM,'NM');
                
            case 'run_NM'
            
                load(pNM);
                [ID, NUM] = nk_CountUniqueOccurences(NM.label);
                minF_CV2 = min(NUM);
                minF_CV1 = min(NUM)-1;
                NM.TrainParam.RAND.OuterFold = minF_CV2;
                NM.TrainParam.RAND.OuterPerm = 1;
                NM.TrainParam.RAND.InnerFold = minF_CV1;
                NM.TrainParam.RAND.InnerPerm = 1;
                GridAct = true(NM.TrainParam.RAND.OuterPerm, NM.TrainParam.RAND.OuterFold );
                NM.analysis{analdim} = nk_TrainML(NM, analdim, 1, GridAct);
                NM = nk_OOCVPrep2(NM, analdim, 1, GridAct, O);
                hu = findobj('Tag','OOCV');
                if ~isempty(hu), delete(hu); end
                save(pNM,'NM');
                
            case 'get_results_NM'
                 fprintf('\nLoading %s',pNM);
                 load(pNM);
                 TP_CV(j,i) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.TP;
                 FP_CV(j,i) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.FP;
                 TN_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.TN;
                 FN_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.FN;
                 Sens_CV(j,i) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.sens;
                 Spec_CV(j,i) = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.spec;
                 BAC_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.BAC;
                 PPV_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.PPV;
                 NPV_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.NPV;
                 PSI_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.PSI;
                 pLR_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.pLR;
                 DOR_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.DOR;
                 NND_CV(j,i)  = NM.analysis{analdim}.GDdims{1}.BinClass{1}.contigency.NND;
                 TP_OOCV(j,i) = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.TP;
                 FP_OOCV(j,i) = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.FP;
                 TN_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.TN;
                 FN_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.FN;
                 Sens_OOCV(j,i) = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.sens;
                 Spec_OOCV(j,i) = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.spec;
                 BAC_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.BAC;
                 PPV_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.PPV;
                 NPV_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.NPV;
                 PSI_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.PSI;
                 pLR_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.pLR;
                 DOR_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.DOR;
                 NND_OOCV(j,i)  = NM.analysis{analdim}.OOCV{1}.BinResults{1}.contingency.NND;
        end
        
    end
    
end
if strcmp(act,'get_results_NM')
    save Perm_Results.mat TP_CV FP_CV TN_CV FN_CV Sens_CV Spec_CV BAC_CV PPV_CV NPV_CV PSI_CV pLR_CV DOR_CV NND_CV ...
        TP_OOCV FP_OOCV TN_OOCV FN_OOCV Sens_OOCV Spec_OOCV BAC_OOCV PPV_OOCV NPV_OOCV PSI_OOCV pLR_OOCV DOR_OOCV NND_OOCV
    H = {'TP','FP','TN','FN','Sens','Spec','BAC','PPV','NPV' ,'PSI','pLR','DOR','NND'};
    DOR_CV(isinf(DOR_CV)) = NaN; 
    pLR_CV(isinf(pLR_CV)) = NaN;
    NND_CV(isinf(NND_CV)) = NaN;
    A = nm_nanmean(TP_CV);            B = nm_nanstd(TP_CV);
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
    A = A';                 B = B';
   
    DOR_OOCV(isinf(DOR_OOCV)) = NaN; 
    pLR_OOCV(isinf(pLR_OOCV)) = NaN;
    NND_OOCV(isinf(NND_OOCV)) = NaN;
    C = nm_nanmean(TP_OOCV);            D = nm_nanstd(TP_OOCV);
    C = [C; nm_nanmean(FP_OOCV)];       D = [D; nm_nanstd(FP_OOCV)];
    C = [C; nm_nanmean(TN_OOCV)];       D = [D; nm_nanstd(TN_OOCV)];
    C = [C; nm_nanmean(FN_OOCV)];       D = [D; nm_nanstd(FN_OOCV)];
    C = [C; nm_nanmean(Sens_OOCV)];     D = [D; nm_nanstd(Sens_OOCV)];
    C = [C; nm_nanmean(Spec_OOCV)];     D = [D; nm_nanstd(Spec_OOCV)]; 
    C = [C; nm_nanmean(BAC_OOCV)];      D = [D; nm_nanstd(BAC_OOCV)];
    C = [C; nm_nanmean(PPV_OOCV)];      D = [D; nm_nanstd(PPV_OOCV)];
    C = [C; nm_nanmean(NPV_OOCV)];      D = [D; nm_nanstd(NPV_OOCV)];
    C = [C; nm_nanmean(PSI_OOCV)];      D = [D; nm_nanstd(PSI_OOCV)];
    C = [C; nm_nanmean(pLR_OOCV)];      D = [D; nm_nanstd(pLR_OOCV)];
    C = [C; nm_nanmean(DOR_OOCV)];      D = [D; nm_nanstd(DOR_OOCV)];
    C = [C; nm_nanmean(NND_OOCV)];      D = [D; nm_nanstd(NND_OOCV)];
    C = C';                 D = D';
    
    Mean_CV_tbl = array2table(A,'VariableNames',H);
    SD_CV_tbl = array2table(B,'VariableNames',H);
    Mean_OOCV_tbl = array2table(C,'VariableNames',H);
    SD_OOCV_tbl = array2table(D,'VariableNames',H);
    writetable(Mean_CV_tbl,'Perm_Results.xls','Sheet','Mean_CV');
    writetable(SD_CV_tbl,'Perm_Results.xls','Sheet','SD_CV');
    writetable(Mean_OOCV_tbl,'Perm_Results.xls','Sheet','Mean_OOCV');
    writetable(SD_OOCV_tbl,'Perm_Results.xls','Sheet','SD_OOCV');
end
clear NM