function visdata = nk_VisModels(inp, id, GridAct, batchflag)
% =========================================================================
% visdata = nk_VisModels(inp, strout, id, GridAct, batchflag)
% =========================================================================
% Visualisation module of NeuroMiner
% This core NeuroMiner function analyzes the relevance of the elements in 
% the predictive patterns of the models and optionally computed model 
% significance using a permutation-based approach
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2020

global SVM SAV RFE MODEFL CV VERBOSE FUSION MULTILABEL EVALFUNC CVPOS

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
visdata         = [];                               % Initialize with empty output
switch inp.analmode
    case 0
        ovrwrt  = inp.ovrwrt;                       % overwrite existing data files
    case 1
        vismat  = inp.vismat;                       % Visualization datamat
end
saveparam       = inp.saveparam;                    % Save parameters to disk
CV1op           = inp.CV1;                          % Operate preprocessing module at the CV1 level to save RAM
loadparam       = inp.loadparam;                    % Use existing parameters loaded from disk
strout          = inp.varstr;                       % Suffix for filename indicating modality assignment of file
nclass          = inp.nclass;                       % Number of binary comparisons
analysis        = inp.analysis;                     % GDanalysis structure to be used
varind          = inp.tF;                           % Current modality index
ol              = 0;                                % Counter for computing mean vectors
ll              = 1;                                % Counter for looping through CV2 arrays
lx              = size(inp.labels,1);               % number of cases
[ix, jx]        = size(CV.TrainInd);                % No. of Perms / Folds at CV2 level
ind0            = false(1,ix*jx);
algostr         = GetMLType(SVM);                   % Algorithm descriptions
FullPartFlag    = RFE.ClassRetrain;                 % Has the user activated full CV1 partition retraining ?
nM              = numel(inp.tF);                    % Number of modalities with independent preprocessing
decompfl        = false(1,nM);                      % flag for factorization methods during preprocessing
permfl          = false(1,nM);                      % flag for permutation mode 
sigfl           = false(1,nM);                      % flag for significance estimation mode
nperms          = ones(1,nM);                       % number of permuations per modality
pmode           = ones(1,nM);                       % Permutation mode
if strcmp(MODEFL,'classification') && nclass > 1
    ngroups     = numel(unique(inp.labels(~isnan(inp.labels)))); % Number of classes in the label
else
    ngroups     = 1;
end
% Loop through modalities (in early fusion: nM = 1)
Dall = 0;

% Check whether you have to run label imputation
IMPUTE.flag = false;
if iscell(inp.PREPROC), iPREPROC = inp.PREPROC{1}; else iPREPROC = inp.PREPROC; end    
if isfield(iPREPROC,'LABELMOD') && isfield(iPREPROC.LABELMOD,'LABELIMPUTE'); 
    IMPUTE = iPREPROC.LABELMOD.LABELIMPUTE; 
    IMPUTE.flag = true; 
end
linsvmfl = determine_linsvm_flag(SVM);

% Always set to binary preprocessing (unless true multi-group learners have
% been intergrated in NM)
BINMOD = iPREPROC.BINMOD; 
clc

inp.id = id;
CVPOS.fFull = FullPartFlag;

for i = 1 : nM
    
    % Dimensionality of current modality
    D = getD(FUSION.flag, inp, i);
    
    % Dimensionality of the (concatenated feature space)
    Dall = Dall + D;
    
    % Activate preprocessing params of current modality
    switch FUSION.flag
        case 2
            iPREPROC = inp.PREPROC{i}; 
            iVis = inp.VIS{i};
        otherwise
            iPREPROC = inp.PREPROC;
            iVis = inp.VIS;
    end
    
    % Determine if factorization methods are involved in current
    % preprocessing chain
    if isfield(iPREPROC,'ACTPARAM')
        for zu = numel(iPREPROC.ACTPARAM) : -1 : 1
            if strcmp(iPREPROC.ACTPARAM{zu}.cmd,'reducedim') ...
                || strcmp(iPREPROC.ACTPARAM{zu}.cmd,'remvarcomp') 
                decompfl(i) = true; 
            end
        end
    end
    
    % Check whether permutation mode is activated in the visualisation setup of current modality 
    if isfield(iVis,'PERM') && iVis.PERM.flag
        if ~isfield(iVis.PERM,'mode')
            pmode(i) = 1; 
        else
            pmode(i) = iVis.PERM.mode;
        end
        permfl(i)      = true; 
        nperms(i)      = iVis.PERM.nperms;
        compfun        = nk_ReturnEvalOperator(SVM.GridParam);
        
        % Generate or load permutation matrix from file.
        cprintf('black*','\nPERMUTATION MODE ENABLED.')
        if pmode(i)
            permfile = fullfile(inp.rootdir,[SAV.matname '_VISpermmat_ID' id '.mat']);
            if ~exist(permfile,'file')
                fprintf('\nCreating parent permutation matrix with %g perms', nperms(1))
                indpermA = nk_VisXPermHelper('genpermlabel', size(inp.labels,1), nperms(1));
                save(permfile,'indpermA');
            else
                fprintf('\nLoading parent permutation matrix from file: \n%s', permfile);
                load(permfile,'indpermA');
                if size(indpermA,2)~= nperms(1)
                    error('Number of permutations in matrix (%g) does not match number of expected permutations (%g)', size(indpermA,2), nperms(1))
                end
                if size(indpermA,1) ~= size(inp.labels,1) 
                    error('Number of cases in permutation matrix (%g) does not match number of cases in NM file (%g)', size(indpermA,1), size(inp.labels,1))
                end
            end
        end
        % Check whether reconstruction of only those components that are significant is
        % activated by the user.
        if decompfl(i) && isfield(iVis.PERM,'sigflag') && iVis.PERM.sigflag,
            sigfl(i) = true;
            sigflFDR = false;
        end
    end
    if isfield(iVis,'use_template') && iVis.use_template,
        fprintf('\nPerform template processing!')
        templateflag(i) = true;
    end        
end

% Initialize CV2 data containers
I = nk_VisXHelper('init', nM, nclass,  decompfl, permfl, ix, jx);
if any(permfl)
    I.VCV2MPERM_S = cell(lx,nclass,nperms(1)); 
    I.VCV2MORIG_S = cell(lx,nclass); 
    if inp.multiflag
        I.VCV2MORIG_S_MULTI = nan(lx,ix);
        I.VCV2MPERM_S_MULTI = nan(lx,ix, nperms(1));
        I.VCV2MPERM_MULTI   = zeros(1,ix*jx);
    end
end

% Obtain feature labels for the selected modalities
featnames = []; if isfield(inp,'featnames'), featnames = inp.featnames; end

% Multi-Group processing ?
multiflag = false; if isfield(inp,'multiflag') && ~isempty(inp.multiflag), multiflag = inp.multiflag; end
if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end
multlabelstr = '';  if MULTILABEL.flag, multlabelstr = sprintf('t%g',inp.curlabel); end

 % Do we have to scale the labels?
[ inp ] = nk_ApplyLabelTransform( inp.PREPROC, MODEFL, inp );

% %%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds
        
        fprintf('\n---------------------------------------------------------------------------------------------')
        if ~GridAct(f,d), 
            ll=ll+1;
            fprintf('\nSkipping CV2 [%g,%g] (user-defined).',f,d)
            continue 
        end;
        
        [iy, jy] = size(CV.cvin{f,d}.TrainInd); % No. of Perms / Folds at CV1 level
        CVPOS.CV2p = f;
        CVPOS.CV2f = d;
        operm = f; ofold = d;
        oVISpath = nk_GenerateNMFilePath(inp.rootdir, SAV.matname, inp.datatype, multlabelstr, strout, id, operm, ofold);
        OptModelPath = nk_GenerateNMFilePath( inp.saveoptdir, SAV.matname, 'OptModel', [], inp.varstr, inp.id, operm, ofold);
    
        switch inp.analmode 
            
            case 0
        
                %%%%%%%%%%%%%%%%%%%%%%%%% USE PRECOMPUTED ? %%%%%%%%%%%%%%%%%%%%%%%
                
                if exist(oVISpath,'file') && ~ovrwrt && ~batchflag
                    
                    [~, onam] = fileparts(oVISpath);
                    fprintf('\nVISdatamat found for CV2 [%g,%g]:',f,d)
                    fprintf('\nLoading: %s',onam)
                    [I, I1, filefound] = nk_VisXHelper('accum', nM, nclass, decompfl, permfl, ix, jx, I, inp, ll, nperms, oVISpath);
                    if filefound
                        if isfield(I1,'PCV1SUM'), PCV1SUMflag=true; else PCV1SUMflag = false; end
                        if any(permfl)
                            for h=1:nclass
                                switch MODEFL
                                    case 'classification'
                                        TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                                    case 'regression'
                                        TsInd2 = CV.TestInd{f,d};
                                end
                                I.VCV2MORIG_S(TsInd2,h) = cellmat_mergecols(I.VCV2MORIG_S(TsInd2,h), num2cell(I1.DS{h},2));
                                for perms = 1:nperms(1), I.VCV2MPERM_S(TsInd2,h,perms) = cellmat_mergecols(I.VCV2MPERM_S(TsInd2,h,perms), num2cell(I1.DS_perm{h}(:,:,perms),2)); end
                            end
                        end
                        ll=ll+1;
                        ind0(ll) = true;
                        ol=ol+1; continue
                    end
                    
                elseif exist(oVISpath,'file') && batchflag
                    
                    % in batch mode we do not compute statistics across the
                    % CV2 partitions
                    [~, onam] = fileparts(oVISpath);
                    fprintf('\nVISdatamat found for CV2 [%g,%g]:\n%s',f,d,onam)
                    fprintf('\nBatch mode detected. Continue.')
                    ll=ll+1;
                    continue

                end

                %%%%%%%%% GET PREPROCESSING PARAMETERS FOR CUR. CV2 PART. %%%%%%%%%
                % First generate parameter array for preprocessing based on
                % the trained base learners in the ensemble. This saves
                % computational resources because we are not going to preprocess the 
                % data with all possible parameter combinations specified by the user, 
                % but only with those chosen by the NM training process.
                
                inp.f = f; inp.d = d; inp.ll = ll;  
                % Parameter flag structure for preprocessing
                paramfl = struct('use_exist',   loadparam, ...
                                 'found',       false, ... 
                                 'write',       true, ... % has to be set to true otherwise no params will be returned from the preproc module
                                 'CV1op',       CV1op, ...
                                 'multiflag',   multiflag);
                
                % Apply prerpocessing on the entire data and use these
                % parameters to adjust for arbitrary PCA rotations through 
                % the Procrustes transform 
                if exist('templateflag','var') && any(templateflag), paramfl.templateflag = true; end
                               
                % Compute params
                inp.loadGD = true;
                if isfield(inp,'CV1') && inp.CV1 == 1, inp.smoothonly = true; end
                [ inp, contfl, analysis, mapY, GD, MD, Param, paramfl ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl);
                inp.loadGD = false;
                
                if contfl, continue; end
                    
                % Prepare containers & initialize matrices for CV1-level
                % weight vector relevance metrics
                % =========================================================
                ol                                = ol+1;
                [~, I1] = nk_VisXHelper('initI1', nM, nclass, decompfl, permfl);
                GDFEAT                         = GD.FEAT; 
                GDVI                           = GD.VI; 
                if inp.stacking
                    if strcmp(SVM.prog,'SEQOPT')
                        mChnl = ones(1,numel(GD.nM_cnt));
                    else
                        mChnl = GD.nM_cnt;
                    end
                end
                clear GD
                
                % Try to load models from disk if user
                % chose to do this
                fndMD = false; 
                if loadparam && isfield(inp,'optmodelmat') && exist(inp.optmodelmat{operm,ofold},'file')
                    fprintf('\nLoading OptModel: %s', inp.optmodelmat{operm,ofold});
                    load(inp.optmodelmat{operm,ofold},'MD'); fndMD = true; 
                end
                if ~fndMD, MD = cell(nclass,1); end
                
                % ---------------------------------------------------------
                if ~VERBOSE,fprintf('\n\nComputing visualizations for CV2 [ %g, %g ] ',f,d), end
                
                %% Initialize containers for analysis
                if permfl, 
                    I1.TS           = cell(nclass,1);
                    I1.DS           = cell(nclass,1);
                    I1.TS_perm      = cell(nclass,1);
                    I1.DS_perm      = cell(nclass,1);
                    if inp.multiflag
                        I1.mTS      = cell(nclass,1);
                        I1.mDS      = cell(nclass,1);
                        I1.mDS_perm = cell(nclass,1);
                        I1.mTS_perm = cell(nclass,1);
                    end
                end
                
                for h=1:nclass % Loop through binary comparisons
                    
                    switch MODEFL
                        case 'classification'
                            TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                            if inp.multiflag, TsIndM = CV.TestInd{f,d}; end
                        case 'regression'
                            TsInd2 = CV.TestInd{f,d};
                    end
                    
                     %% Step 1: Get optimal model parameters
                    % Retrieve optimal parameters from precomputed analysis structure
                    % Differentiate according to binary or multi-group mode
                    [~, Pspos, nP] = nk_GetModelParams2(analysis, multiflag, ll, h);
                    
                    % Allocate memory to store CV1 ensemble patterns
                    ill = getModelNumDim(h,iy,jy,nP,Pspos,GDFEAT);
                    fprintf('\nPredictor #%g: Need to evaluate %g models in this CV2 partition', h, ill);
                    
                    if permfl, 
                        I1.VCV1WPERM{h} = cell(ill,1); 
                        I1.VCV1MPERM{h} = nan(ill,1); 
                    end
                    
                    for n=1:nM
                        
                        % Retrieve dimensionality of target space
                        D                          = getD(FUSION.flag, inp, n);
                        % Setup container for weight storage
                        I1.VCV1{h,n}               = zeros(D, ill,'single'); 
                        I1.VCV1SUM{h,n}            = zeros(D, 1,'single');
                        I1.VCV1SQ{h,n}             = I1.VCV1SUM{n};
                        I1.VCV1MEAN{h,n}           = I1.VCV1SUM{n};
                        I1.VCV1STD{h,n}            = I1.VCV1SUM{n};
                        
                        % Prepare for analysis without factorization
                        if any(~decompfl), 
                            I1.PCV1SUM{h, n}                    = zeros(D, 1,'single'); 
                            I1.VCV1PEARSON{h, n}                = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1SPEARMAN{h, n}               = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1PEARSON_UNCORR_PVAL{h, n}    = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}   = nan(D, iy*jy*nP,'single');
                            I1.VCV1PEARSON_FDR_PVAL{h, n}       = nan(D, iy*jy*nP,'single'); 
                            I1.VCV1SPEARMAN_FDR_PVAL{h, n}      = nan(D, iy*jy*nP,'single'); 
                        end
                        
                        % Prepare for permutation analysis
                        if permfl
                            nTs                 = size(TsInd2,1);
                            I1.TS{h}            = nan(nTs, ill);
                            I1.DS{h}            = nan(nTs, ill);
                            I1.DS_perm{h}       = nan(nTs, ill, nperms(1));
                            I1.TS_perm{h}       = nan(nTs, ill, nperms(1));
                            if inp.multiflag
                                nTsM = size(CV.TestInd{f,d},1);
                                I1.mTS{h}       = nan(nTsM, ill);
                                I1.mDS{h}       = nan(nTsM, ill);
                                I1.mDS_perm{h}  = nan(nTsM, ill, nperms(1));
                                I1.mTS_perm{h}  = nan(nTsM, ill, nperms(1));
                            end
                            I1.VCV1PERM{h, n}   = nan(D, ill,'single');
                            I1.VCV1PERM_FDR{h, n} = nan(D, ill,'single');
                            I1.VCV1ZSCORE{h, n} = nan(D, ill,'single');
                        end
                    end
                    
                    % Initialize model container 
                    if ~fndMD , MD{h} = cell(nP,1); for m = 1 : nP, MD{h}{m} = cell(iy,jy); end; end;
                end
                
                il = ones(nclass,1); kil=ones(nclass,1);
                
                %% Perform analyses   
                for k=1:iy % CV1 permutations

                    for l=1:jy % CV1 folds
                        
                        CVPOS.CV1p = k;
                        CVPOS.CV1f = l;
                        
                        if isfield(inp,'CV1') && inp.CV1 == 1
                            inp.CV1p = [k,k]; inp.CV1f = [l,l];
                            fprintf('\nPreprocessing data at selected parameter combination(s) ');
                            [ ~, ~, analysis, mapY, ~, ~, Param, paramfl ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl);
                        end
                        
                        for h=1:nclass % Loop through binary comparisons
                    
                            if nclass > 1, fprintf('\n');cprintf('blue*', '*** %s #%g *** ',algostr, h); end

                            switch MODEFL
                                case 'classification'
                                    TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                                    if inp.multiflag, TsIndM = CV.TestInd{f,d}; end
                                case 'regression'
                                    TsInd2 = CV.TestInd{f,d};
                            end
                            
                            TsInd = mapY.TsInd{h}; 
                            
                            %% Step 1: Get optimal model parameters
                            % Retrieve optimal parameters from precomputed analysis structure
                            % Differentiate according to binary or multi-group mode
                            [Ps, Pspos, nP, Pdesc] = nk_GetModelParams2(analysis, multiflag, ll, h);
                            
                            for m = 1 : nP % parameter combinations

                                if nP>1, fprintf('\n');cprintf(rgb('RoyalBlue'),'Extracing model parameters at parameter node #%g/%g ', m, nP); end
                                % Prepare learning params
                                cPs = Ps(m,:); sPs = nk_PrepMLParams(Ps, Pdesc, m);

                                % -----------------------------------------------------
                                % Construct pattern for every base learnern in
                                % current CV1 [k,l] partition:
                                %% CV1 LOOP
                                P_str = nk_DefineMLParamStr(cPs, analysis.Model.ParamDesc, h);
                                
                                Fkl = GDFEAT{Pspos(m)}{k,l,h}; 
                                 
                                % Determine number of features in mask and
                                % convert feature mask to logical index
                                ul=size(Fkl,2);
                                if ~islogical(Fkl),F = Fkl ~= 0; else F = Fkl; end

                                if VERBOSE
                                    fprintf('\n');cprintf('*black',['Constructing predictive pattern(s) in CV2 [ %2g ,%2g ], ' ...
                                    'CV1 [ %2g ,%2g, Predictor #%g/%g ]: %g model(s), %s ML params [ %s ]. '], f, d, k, l, h, nclass, ul, algostr, P_str); 
                                else
                                    fprintf('\n');cprintf('*black','Visualizing: CV2 [ %2g, %2g ], CV1 [ %2g, %2g, P: #%g/%g ]: %g model(s) ',f, d, k, l, h, nclass, ul) ;
                                end

                                CVInd   = mapY.CVInd{k,l}{h};
                                TrInd   = mapY.TrInd{k,l}{h};
                                                        
                                % Set the pointer to the correct mapY shelf
                                for n=1:numel(paramfl)
                                    pnt = 1;
                                    if ~BINMOD
                                         if isfield(paramfl{n},'PREPROC') && ...
                                           isfield(paramfl{n},'PXfull') && ...
                                           ~isempty(paramfl{n}.P{1})
                                            pnt = m;
                                            break   
                                        end
                                    else
                                        if isfield(paramfl{n},'PREPROC') && ...
                                           isfield(paramfl{n},'PXfull') && ...
                                           ~isempty(paramfl{n}.P{h})
                                            %Actual parameter node
                                            pnt = m; 
                                            break   
                                        end
                                    end
                                end
                                
                                %%%%%% RECOMPUTE ORIGINAL MODEL %%%%%%
                                % get CV1 training and test data
                                if BINMOD, hix = h; else, hix =1; end
                                [ modelTr , modelCV, modelTs] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, [], Param{1}(k,l,hix), pnt); 
                                switch FUSION.flag
                                    case 2
                                        ParamX = cell(n,1);
                                        for n=1:nM
                                            [~,~,~,~, ParamX{n} ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, [], Param{n}(k,l,hix), pnt); 
                                        end
                                    otherwise
                                         [~,~,~,~, ParamX ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, [], Param{1}(k,l,hix), pnt); 
                                end      

                                % Get label info
                                modelTrL = mapY.TrL{k,l}{h};
                                modelCVL = mapY.CVL{k,l}{h};
                                
                                % Impute labels if needed
                                [modelTrL, modelTr, TrInd] = nk_LabelImputer( modelTrL, modelTr, TrInd, sPs, IMPUTE);
                                
                                % Concatenate Training and CV data if needed
                                if FullPartFlag, 
                                    modelTr = [modelTr; modelCV ]; 
                                    modelTrL = [modelTrL; modelCVL]; 
                                    TrInd = [TrInd; CVInd]; 
                                end
                                modelTr = modelTr(TrInd,:);
                                
                                % Prepare permutation operations
                                if any(permfl),
                                    indperm = []; indpermfeat = [];
                                    switch MODEFL
                                        case 'classification'
                                            pTrInd = CV.TrainInd{f,d}(CV.class{f,d}{h}.TrainInd{k,l}); 
                                            pCVInd = CV.TrainInd{f,d}(CV.class{f,d}{h}.TestInd{k,l}); 
                                        otherwise
                                            pTrInd = CV.TrainInd{f,d}(CV.cvin{f,d}.TrainInd{k,l}); 
                                            pCVInd = CV.TrainInd{f,d}(CV.cvin{f,d}.TestInd{k,l});
                                    end
                                    if size(pTrInd,2)>1,pTrInd=pTrInd'; end
                                    if size(pCVInd,2)>1,pCVInd=pCVInd'; end
                                    if FullPartFlag, pTrInd = [pTrInd; pCVInd];end
                                    if pmode(1), indperm = indpermA(pTrInd,:); end
                                    if inp.multiflag ==1,
                                        modelTsm = modelTs;
                                        modelTsmL = inp.labels(TsIndM);
                                    end
                                    modelTs = modelTs(TsInd,:); 
                                    modelTsL = mapY.TsL{h};
                                else
                                    modelTs = []; modelTsL = [];
                                end
                                
                                if ~fndMD, MD{h}{m}{k,l} = cell(ul,1); end
                                
                                % Loop through base learners' feature masks
                                for u=1:ul
                                     
                                    % Extract features according to mask
                                    try
                                        Ymodel = nk_ExtractFeatures(modelTr, F, [], u);
                                    catch
                                        error('Dimensionality mismatch between training data matrix and feature selection mask. Check your settings')
                                    end
                                    Find = F(:,u);
                                    %If permutation mode expects feature
                                    % permutation prepare for this here:
                                    if any(permfl) && pmode(1) > 1, indpermfeat = nk_VisXPermHelper('genpermfeats', sum(F(:,u)), nperms(1), modelTrL); end
                                    
                                    if ~isempty(GDVI{Pspos(m)})
                                        Vind = GDVI{Pspos(m)}{k,l,h}(:,u);
                                    else
                                        Vind = true(size(Find,1),1);
                                    end
                                    
                                    % Model computation
                                    if isfield(mapY,'VI')
                                        if ~fndMD, [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(Ymodel, modelTrL, sPs, 1, mapY.VI{k,l}{hix}{u}); end
                                    else
                                        if ~fndMD, [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(Ymodel, modelTrL, sPs, 1); end
                                    end
                                    
                                    if inp.stacking
                                        vec_mj = [];
                                        for mj = 1:inp.nD
                                            vec_mj = [vec_mj; mj*ones(mChnl(mj),1)];
                                        end
                                    end
                                    
                                    if ~any(permfl)
                                        % Compute original weight map in input space
                                        [Tx, Psel, Rx, SRx, ~, PAx] = nk_VisXWeight2(inp, MD{h}{m}{k,l}{u}, Ymodel, modelTrL, varind, ParamX, Find, Vind, decompfl);
                                    else
                                        fprintf('\n\t%3g | OptModel =>', u);
                                         % Build permutation index array if
                                        % needed and compute further variables
                                        % needed for perm-based stats
                                        if inp.multiflag
                                            [~, I1.mTS{h}(:,il(h)), I1.mDS{h}(:,il(h))] = nk_GetTestPerf(modelTsm, modelTsmL, Find, MD{h}{m}{k,l}{u}, modelTr); 
                                        end
                                        [perf_orig, I1.TS{h}(:,il(h)), I1.DS{h}(:,il(h))] = nk_GetTestPerf(modelTs, modelTsL, Find, MD{h}{m}{k,l}{u}, modelTr); 
                                        fprintf(' %1.2f',perf_orig)
                                        % if sigfl = true
                                        % Determine significant components
                                        % through non-parametric permutation
                                        if sigfl
                                            fprintf(' | Significant pattern components (%g perms):\t',nperms(1));
                                            % Original weight vector
                                            [~,~,~,~,Vx] = nk_VisXWeight2(inp, MD{h}{m}{k,l}{u}, Ymodel, modelTrL, varind, ParamX, Find, Vind, decompfl, false);
                                            ipVx = Vx >= 0; inVx = Vx < 0;
                                            Vx_perm = zeros( size(Vx,1), nperms(1) );
                                            MD_perm = cell(nperms(1),1);
                                            for perms = 1:nperms(1)
                                                fprintf('+');
                                                % Train permuted model
                                                [L_perm, Ymodel_perm]       = nk_VisXPermY(Ymodel, inp.labels, pTrInd, pmode(1), indperm, indpermfeat, perms);
                                                if isfield(mapY,'VI')
                                                    [~, MD_perm{perms}]     = nk_GetParam2(Ymodel_perm, L_perm, sPs, 1, mapY.VI{k,l}{hix}{u});
                                                else
                                                    [~, MD_perm{perms}]     = nk_GetParam2(Ymodel_perm, L_perm, sPs, 1);
                                                end
                                                % Compute permuted weight vector in feature space
                                                [~,~,~,~,Vx_perm(:, perms)] = nk_VisXWeight2(inp, MD_perm{perms}, Ymodel_perm, L_perm, varind, ParamX, Find, Vind, decompfl, false);
                                            end
                                            % Determine significance by comparing observed weight vector
                                            % against null distribution of permuted models' weight vectors
                                            VxV = zeros(size(Vx));
                                            VxV(ipVx) = (sum(bsxfun(@ge,Vx_perm(ipVx,:),Vx(ipVx)),2)/nperms(1))/ul;
                                            VxV(inVx) = (sum(bsxfun(@le,Vx_perm(inVx,:),Vx(inVx)),2)/nperms(1))/ul;
                                            I1.VCV1WPERM{il(h),h} = VxV; 
                                            clear VxV
                                            if sigflFDR
                                                Fadd   = fdr_bh(I1.VCV1WPERM{il(h),h},0.5,'dep');
                                                FDRstr = '(FDR) ';
                                            else
                                                Fadd   = (I1.VCV1WPERM{il(h),h} <= 0.5);
                                                FDRstr = '(uncorr) ';
                                            end
                                            if ~sum(Fadd)
                                                [minP, Fadd] = min(I1.VCV1WPERM{il(h),h});
                                                fprintf('\tNo component significant at alpha %s = 0.5 => relaxing to max P = %g\n\t\t\t\t\t\t\t\t',FDRstr, minP );
                                            else
                                                fprintf('\t%g / %g components significant at alpha %s = 0.5\n\t\t\t\t\t\t\t\t',sum(Fadd), numel(Fadd),FDRstr);
                                            end
                                        else
                                            Fadd = true(size(F,1),1);
                                        end
                                        
                                        % Compute original weight map in input space
                                        [Tx, Psel, Rx, SRx, ~, PAx ] = nk_VisXWeight2(inp, MD{h}{m}{k,l}{u}, Ymodel, modelTrL, varind, ParamX, Find, Vind, decompfl, [], Fadd);
                                        fprintf(' | Permuting:\t');
                                        Tx_perm = cell(1,nM); Px_perm = zeros(1,nperms(1));
                                        for n=1:nM, Tx_perm{n} = zeros(size(Tx{n},1),nperms(n)); end
                                        for perms = 1:nperms(1)
                                            if ~sigfl, 
                                                % Train permuted model
                                                [ L_perm, Ymodel_perm ] = nk_VisXPermY(Ymodel, inp.labels, pTrInd, pmode(1), indperm, indpermfeat, perms);
                                                  if isfield(mapY,'VI')
                                                      [~, MDs] = nk_GetParam2(Ymodel_perm, L_perm, sPs, 1, mapY.VI{k,l}{hix}{u} );
                                                  else
                                                      [~, MDs] = nk_GetParam2(Ymodel_perm, L_perm, sPs, 1 );
                                                  end
                                            else
                                                % Retrieve trained permuted model
                                                MDs = MD_perm{perms};
                                            end
                                            % Compute permuted model test performance
                                            [perf_perm, I1.TS_perm{h}(:,il(h),perms), I1.DS_perm{h}(:,il(h),perms)] = nk_GetTestPerf(modelTs, modelTsL, Find, MDs, modelTr);
                                            if inp.multiflag
                                                [~, I1.mTS_perm{h}(:,il(h),perms), I1.mDS_perm{h}(:,il(h),perms)] = nk_GetTestPerf(modelTsm, modelTsmL, Find, MDs, modelTr); 
                                            end
                                            % Compare against original model performance
                                            if feval(compfun, perf_perm, perf_orig),
                                                fprintf('.'); 
                                                Px_perm(perms) = Px_perm(perms) + 1;
                                            end
                                            % Compute permuted weight map in input space
                                            TXperms = nk_VisXWeight2(inp, MDs, Ymodel_perm, L_perm, varind, ParamX, Find, Vind, decompfl, [], Fadd);
                                            for n=1:nM, Tx_perm{n}(:,perms) = TXperms{n}; end
                                        end
                                        
                                        % Model significance
                                        I1.VCV1MPERM{h}(il(h)) = (sum(Px_perm) / nperms(1));
                                        
                                        % Print significance to screen
                                        cprintf('black*', ' P=%1.3f ', I1.VCV1MPERM{h}(il(h)));
                                        
                                        % Loop through modalities
                                        for n=1:nM
                                            % Define index vector to
                                            % original space of modality
                                            Fpind = any(Tx_perm{n},2)'; 
                                            
                                            % Now compute the P value vector:
                                            TxV = zeros(size(Tx{n}));
                                            ipTx = Tx{n} >= 0; inTx = Tx{n} < 0;
                                            TxV(ipTx) = sum(bsxfun(@ge,Tx_perm{n}(ipTx,:),Tx{n}(ipTx)),2)/nperms(1);
                                            TxV(inTx) = sum(bsxfun(@le,Tx_perm{n}(inTx,:),Tx{n}(inTx)),2)/nperms(1);
                                            Pvals = TxV(Fpind);
                                            
                                            % ... and the Z score vector:
                                            Zvals = bsxfun(@rdivide, Tx{n} - nm_nanmean(Tx_perm{n},2), nm_nanstd(Tx_perm{n},2));
                                            Zvals = Zvals(Fpind);
                                            
                                            % with stacking there is a
                                            % different way to compute
                                            % significance, as we have to
                                            % find to which child predictor a
                                            % feature belongs.
                                            if inp.stacking
                                                for mj = 1:inp.nD
                                                    mjPvals = mean(Pvals(vec_mj(Fpind) == mj));
                                                    mjC = mean(Zvals(vec_mj(Fpind) == mj));
                                                    I1.VCV1PERM{h,n}(mj,il(h)) = mjPvals;
                                                    [~,~,~,I1.VCV1PERM_FDR{h,n}(mj, il(h))] = fdr_bh(mjPvals,0.05,'pdep');
                                                    I1.VCV1ZSCORE{h,n}(mj,il(h)) = mjC;
                                                end
                                            else
                                                [ ~, ~, ~, badcoords ] = getD(FUSION.flag, inp, n); badcoords = ~badcoords;
                                                I1.VCV1PERM{h,n}(badcoords & Fpind,il(h)) = Pvals;
                                                [~,~,~,I1.VCV1PERM_FDR{h,n}(badcoords & Fpind, il(h))] = fdr_bh(Pvals,0.05,'pdep'); 
                                                I1.VCV1ZSCORE{h,n}(badcoords & Fpind,il(h)) = Zvals;
                                            end
                                            % and show how many uncorrected P
                                            % values are below alpha=0.05
                                            sigcomp = sum(I1.VCV1PERM{h,n}(:,il(h))<=0.05);
                                            sigmin  = min(I1.VCV1PERM{h,n}(:,il(h)));
                                            FDRsigmin = min(I1.VCV1PERM_FDR{h,n}(:,il(h)));
                                            fprintf('; Modality #%g: [ %g features <= 0.05, min P value (uncorr, FDR) = %1.6f, %1.6f ] ', n, sigcomp, sigmin, FDRsigmin);
                                        end
                                    end
                                    
                                    % Some additional computation if
                                    % factorization methods have not been used
                                    % in the preprocessing chain of a modality
                                    for n = 1:nM
                                        Fpind = any(Tx{n},2); 
                                        if inp.stacking
                                            for mj = 1:inp.nD
                                                Imj = vec_mj == mj & Fpind;
                                                I1.VCV1{h,n}(mj,il(h)) = mean(Tx{n}(Imj));
                                                if ~decompfl(n) && u==1,  
                                                    I1.VCV1PEARSON{h, n}(mj,kil(h)) = nm_nanmean(Rx{n}(Imj));
                                                    I1.VCV1SPEARMAN{h, n}(mj,kil(h)) = nm_nanmean(SRx{n}(Imj));
                                                    I1.VCV1PEARSON_UNCORR_PVAL{h, n}(mj,kil(h)) = nm_nanmean(nk_PTfromR(Rx{n}(Imj), size(Ymodel,1), 2));
                                                    I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}(mj,kil(h)) = nm_nanmean(nk_PTfromR(SRx{n}(Imj), size(Ymodel,1), 2));
                                                    if linsvmfl
                                                        I1.VCV1PVAL_ANALYTICAL{h, n}(mj,kil(h)) = nm_nanmean(PAx{n}(Imj));
                                                    end
                                                end
                                            end
                                             %% Comnpute feature selection probabilities %%
                                             if inp.fixedOrder && ~decompfl(n)
                                                 if isempty(I1.PCV1SUM{h, n}), I1.PCV1SUM{h, n} = zeros(inp.nD,1); end
                                                 I1.PCV1SUM{h, n} = I1.PCV1SUM{h, n} + Psel{n}; 
                                             end
                                        else
                                            % Define index vector to
                                            % original space of modality
                                            [ ~, ~, ~, badcoords] = getD(FUSION.flag, inp, n); badcoords = ~badcoords;

                                            % Store results in CV1 container variables                                    
                                            % I1.numCV1parts(h, n) = I1.numCV1parts(h, n) + 1;
                                            I1.VCV1{h,n}(badcoords,il(h)) = Tx{n};

                                            if ~decompfl(n) 
                                                if u==1,  
                                                    %% Compute univariate correlation coefficient for each feature
                                                    I1.VCV1PEARSON{h, n}(badcoords,kil(h)) = Rx{n};
                                                    I1.VCV1SPEARMAN{h, n}(badcoords,kil(h)) = SRx{n};
                                                    I1.VCV1PEARSON_UNCORR_PVAL{h, n}(badcoords,kil(h)) = nk_PTfromR(Rx{n}, size(Ymodel,1), 2);
                                                    I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}(badcoords,kil(h)) = nk_PTfromR(SRx{n}, size(Ymodel,1), 2);
                                                    [~,~,~,I1.VCV1PEARSON_FDR_PVAL{h, n}(badcoords,kil(h))] = fdr_bh(I1.VCV1PEARSON_UNCORR_PVAL{h, n}(badcoords,kil(h)), 0.05, 'pdep');
                                                    [~,~,~,I1.VCV1SPEARMAN_FDR_PVAL{h, n}(badcoords,kil(h))] = fdr_bh(I1.VCV1SPEARMAN_UNCORR_PVAL{h, n}(badcoords,kil(h)), 0.05, 'pdep');
                                                end
                                                %% Comnpute feature selection probabilities %%
                                                if isempty(I1.PCV1SUM{h, n}), I1.PCV1SUM{h, n} = zeros(size(badcoords,2),1); end
                                                I1.PCV1SUM{h, n}(badcoords) = I1.PCV1SUM{h, n}(badcoords) + Psel{n}; 
                                                if linsvmfl
                                                     I1.VCV1PVAL_ANALYTICAL{h, n}(badcoords,kil(h)) = PAx{n};
                                                     try
                                                        [~,~,~,I1.VCV1PVAL_ANALYTICAL_FDR{h, n}(badcoords,kil(h))] = fdr_bh(PAx{n},0.05,'pdep');
                                                     catch
                                                         fprintf('\n')
                                                     end
                                                end
                                            end
                                       end
                                    end
                                    il(h)=il(h)+1;
                                end
                                kil(h)=kil(h)+1;
                                clear Tx Vx Tx_perm Vx_perm tSRx tRx Rx SRx MD_perm
                                %fprintf(' Done.')
                            end
                        end  
                    end
                    if any(permfl)
                        for h=1:nclass
                            % Compute CV2-level model significance
                            modelTsL = mapY.TsL{h};
                            I1.VCV1MORIG_EVALFUNC_CV2{h} = feval(EVALFUNC, modelTsL, nm_nanmedian(I1.DS{h},2)); 
                            I1.VCV1MPERM_CV2{h} = zeros(nperms(1),1); 
                            I1.VCV1MPERM_EVALFUNC_CV2{h} = zeros(nperms(1),1);
                            for perms = 1:nperms(1)
                                I1.VCV1MPERM_EVALFUNC_CV2{h}(perms) = feval(EVALFUNC, modelTsL, nm_nanmedian(I1.DS_perm{h}(:,:,perms),2)); 
                                I1.VCV1MPERM_CV2{h}(perms)          = feval(compfun, I1.VCV1MPERM_EVALFUNC_CV2{h}(perms), I1.VCV1MORIG_EVALFUNC_CV2{h} );
                            end
                        end
                    end
                    clear Tx tmp V Ymodel modelTr modelTrL F Fkl dum 
                end
                I = nk_VisXHelper('accum', nM, nclass, decompfl, permfl, ix, jx, I, inp, ll, nperms(1), I1);  
                fprintf('\nSaving %s', oVISpath); save(oVISpath,'I1','sPs','operm','ofold', '-v7.3');
                if saveparam, fprintf('\nSaving %s', OptModelPath); save(OptModelPath, 'MD', 'ofold','operm', '-v7.3'); end
                if isfield(I1,'PCV1SUM'), PCV1SUMflag=true; else PCV1SUMflag = false; end
                clear Param MD
                
            case 1

                vpth = deblank(vismat{f,d});

                if isempty(vpth) || ~exist(vpth,'file')
                    warning(['No valid VISdata-MAT detected for CV2 partition ' '[' num2str(f) ', ' num2str(d) ']!']);
                else
                    [~,vnam] = fileparts(vpth);
                    ind0(ll) = true; 
                    fprintf('\n\nLoading visualization data:');
                    fprintf('\n%s',vnam);
                    [I, I1] = nk_VisXHelper('accum', nM, nclass, decompfl, permfl, ix, jx, I, inp, ll, nperms(1), vpth);
                    ol=ol+1;
                end
                if isfield(I1,'PCV1SUM'), PCV1SUMflag=true; else, PCV1SUMflag = false; end

        end
        % Assemble observed and permuted model predictions of current CV2
        % partition into overall prediction matrix
        if any(permfl)
            for h=1:nclass
                 fprintf('\nCV2 [%g, %g]: Observed performance [ model #%g ]: %1.2f; P[CV2] = %1.3f', operm, ofold, h, I1.VCV1MORIG_EVALFUNC_CV2{h}, I.VCV2MPERM_CV2(h,ll) ); 
                 if strcmp(MODEFL,'classification')
                    TsInd2 = CV.TestInd{f,d}(CV.classnew{f,d}{h}.ind);
                 else
                    TsInd2 = CV.TestInd{f,d};
                 end
                 I.VCV2MORIG_S(TsInd2,h) = cellmat_mergecols(I.VCV2MORIG_S(TsInd2,h), num2cell(I1.DS{h},2));
                 for perms = 1:nperms(1), I.VCV2MPERM_S(TsInd2,h,perms) = ...
                         cellmat_mergecols(I.VCV2MPERM_S(TsInd2,h,perms), num2cell(I1.DS_perm{h}(:,:,perms),2)); end
            end
            % Perform multi-class permutation testing
            if inp.multiflag
                 if ~isfield(I1,'mDS'), fprintf('\n');error('You requested multi-class significance but the current VISDATAMAT contains only binary classifier data. Rerun visualization with multi-group optimization.'); end
                 fprintf('\nComputing CV2 multi-class model significance\n\tMulti-class perfomance: '); 
                 mDTs = []; mTTs = []; Classes = []; mcolend = 0;
                 TsIndM = CV.TestInd{f,d};
                 for h=1:nclass
                     for il=1:size(I1.mDS{h},2)
                        % Multi-group CV2 array construction for observed
                        % model
                        [mDTs, mTTs, Classes, ~, mcolend] = ...
                            nk_MultiAssemblePredictions( I1.mDS{h}(:,il), I1.mTS{h}(:,il), mDTs, mTTs, Classes, 1, h, mcolend );
                     end
                 end
                 % compute multi-group performance of observed model
                 [ mCV2perf_obs, mCV2pred_obs ] = nk_MultiEnsPerf(mDTs, mTTs, inp.labels(TsIndM, inp.curlabel), Classes, ngroups);
                 I.VCV2MORIG_S_MULTI(TsIndM,f) = mCV2pred_obs;
                 
                 % now assess observed model's significance
                 mPx_perm = zeros(1,nperms(1));
                 fprintf(' %1.2f | Permuting:\t', mCV2perf_obs);
                 for perms = 1:nperms(1)
                     mDTs_perm = []; mTTs_perm = []; Classes_perm = []; mcolend_perm = 0;
                     for h=1:nclass
                         % Multi-group CV2 array construction for permuted
                         % model
                         for il=1:size(I1.mDS{h},2)
                            [mDTs_perm, mTTs_perm, Classes_perm, ~, mcolend_perm] = ...
                                nk_MultiAssemblePredictions( I1.mDS_perm{h}(:,il,perms), I1.mTS_perm{h}(:,il,perms), mDTs_perm, mTTs_perm, Classes_perm, 1, h, mcolend_perm );
                         end
                     end
                     % Compute multi-group performance of permuted
                     % model
                     [ mCV2perf_perm, mCV2pred_perm ] = nk_MultiEnsPerf(mDTs_perm, mTTs_perm, inp.labels(TsIndM, inp.curlabel), Classes_perm, ngroups);
                     I.VCV2MPERM_S_MULTI(TsIndM,f,perms) = mCV2pred_perm;
                     % Compare against original model performance
                     if feval(compfun, mCV2perf_perm, mCV2perf_obs),
                        fprintf('.'); 
                        mPx_perm(perms) = mPx_perm(perms) + 1;
                     end
                 end
                 I.VCV2MPERM_MULTI(ll) = (sum(mPx_perm) / nperms(1));
                 fprintf('P[CV2] = %1.3f', I.VCV2MPERM_MULTI(ll));
            else
                if isfield(I1,'mDS'), fprintf('\nNot computing CV2 multi-class model significance'); end
            end
        end
        if isfield(inp,'issmoothed'), inp.issmoothed = false; end
        ll=ll+1; clear GDFEAT
        clear I1
    end
end


%% PERFORM POST-PROCESSING VISUALIZATION PROCEDURES ON THE ENTIRE DATA 
if ~batchflag
    
    % visdata is a cell array of with nM dimensions, where nM is defined by
    % the number of modalities used in the analysis
    visdata = cell(1,nM);
    
    if any(permfl) 
        
        % mean and SD model significance of binary models at the CV1 level
        I.VCV2MODELP = nm_nanmedian(I.VCV2MPERM,2); 
        I.VCV2MODELP_STD = nm_nanstd(I.VCV2MPERM); 
        
        % Loop through classifiers / regressors and compute hold-out
        % model significance 
        for h=1:nclass
            switch MODEFL
                case 'classification'
                    if numel(CV.class{1,1}{h}.groups) == 2
                        ind1 = inp.labels == CV.class{1,1}{h}.groups(1); f1 = ones(sum(ind1),1);
                        ind2 = inp.labels == CV.class{1,1}{h}.groups(2); f2 = -1*ones(sum(ind2),1);
                        labelh = zeros(numel(inp.labels),1);
                        labelh(ind1) = f1; labelh(ind2) = f2; %labelh(~labelh)=[];
                    else
                        labelh = zeros(size(inp.labels,1),1);
                        ind1 = inp.labels == CV.class{1,1}{h}.groups(1); 
                        labelh(ind1) = 1; labelh(~ind1,h) = -1;
                    end
                case 'regression'
                    labelh = inp.labels;
            end
            indempt     = ~(cellfun(@isempty,I.VCV2MORIG_S) | isnan(cellfun(@sum,I.VCV2MORIG_S)));
            Porig       = cellfun(@nm_nanmedian,I.VCV2MORIG_S(indempt(:,h),h)); 
            Lorig       = labelh(indempt(:,h));
            if inp.targscale, IN.revertflag = true; IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; Porig = nk_PerfScaleObj(Porig, IN); end
            % Observed hold out performance
            I.VCV2MORIG_EVALFUNC_GLOBAL(h) = feval(EVALFUNC, Lorig, Porig);
            fprintf('\nTesting observed %s model performance [ model #%g: %s = %1.2f ] in the hold-out data against %g permutations\n', ...
                                        MODEFL, h, EVALFUNC,  I.VCV2MORIG_EVALFUNC_GLOBAL(h), nperms(1));
            for perms = 1 : nperms(1)
                Pperm                                = cellfun(@nm_nanmedian, I.VCV2MPERM_S(indempt(:,h),h,perms));
                if inp.targscale, Pperm              = nk_PerfScaleObj(Pperm, IN); end
                % Permuted hold-out performances
                I.VCV2MPERM_EVALFUNC_GLOBAL(h,perms) = feval(EVALFUNC, Lorig, Pperm); 
                crt                                  = feval(compfun, I.VCV2MPERM_EVALFUNC_GLOBAL(h,perms), I.VCV2MORIG_EVALFUNC_GLOBAL(h));
                if ~crt, fprintf('*'); else, fprintf('.'); end
                % Store boolean measuring whether permuted hold-out
                % performance was better than the observed performance
                I.VCV2MPERM_GLOBAL(h,perms)          = crt;
            end
        end
        
        % Test additional binary labels
        if ~isempty(inp.extraL)
            I = nk_VisModels_ExtraLabels(I, inp, nperms, compfun);
        end
        
        % [ Eventually in an updated version test subgroup performance
        % here]
        
        % Perform multi-class significance test including entire model and
        % binary dichotomizers (one vs. REST models, irrespective of 
        % one-vs-one or one-vs-all decomposition schemes)
        if inp.multiflag
            
            % mean and SD model significance of multi-class models at the CV1 level
            I.VCV2MODELP = nm_nanmedian(I.VCV2MPERM_MULTI); 
            I.VCV2MODELP_STD = nm_nanstd(I.VCV2MPERM_MULTI); 
            
            I.VCV2MPERM_GLOBAL_MULTI = zeros(1,nperms(1));
            fprintf('\nTesting observed multi-class model performance against permuted models using entire data: %g permutations\n', nperms(1))
            % Convert predictions to class-membership probabilities of
            % observed model
            MultiCV2Prob_orig           = nk_ConvProbabilities(I.VCV2MORIG_S_MULTI, ngroups);
            ind                         = isnan(MultiCV2Prob_orig(:,1)) | isnan(inp.labels); 
            [~, MultiCV2Pred_orig]      = max(MultiCV2Prob_orig,[],2);
            MultiCV2Pred_orig(ind)      = NaN;
            MultiCV2Errs_orig           = nan(size(ind,1),1);
            MultiCV2Errs_orig(ind)      = inp.labels(ind)~= MultiCV2Pred_orig(ind);
            % Compute multi-class confusion matrix of observed model
            MultiCV2ConfMatrix_orig     = nk_ComputeConfMatrix(inp.labels, MultiCV2Pred_orig, ngroups);
            % Observed multi-class performance and one vs. REST
            % dichotomizers are stored in a structure:
            I.VCV2MORIG_GLOBAL_MULTI    = nk_MultiClassAssessConfMatrix(MultiCV2ConfMatrix_orig, inp.labels, MultiCV2Pred_orig, MultiCV2Errs_orig, 'BAC');
            fprintf('\nMulti-class performance: %1.2f | Permuting:\t', I.VCV2MORIG_GLOBAL_MULTI.BAC_Mean);
            I.VCV2PERM_GLOBAL_MULTI     = zeros(1,nperms(1));
            I.VCV2PERM_GLOBAL_MULTI_ONEvsREST = zeros(ngroups,nperms(1));
            for perms = 1 : nperms(1)
                % Convert predictions to class-membership probabilities
                MultiCV2Prob_perm               = nk_ConvProbabilities(I.VCV2MPERM_S_MULTI(:,:,perms), ngroups);
                [~, MultiCV2Pred_perm]          = max(MultiCV2Prob_perm,[],2);
                MultiCV2Pred_perm(ind)          = NaN;
                MultiCV2Errs_perm               = nan(size(ind,1),1);
                MultiCV2Errs_perm(ind)          = inp.labels(ind)~= MultiCV2Pred_perm(ind);
                % Compute multi-class confusion matrix of permuted model
                MultiCV2ConfMatrix_perm         = nk_ComputeConfMatrix(inp.labels, MultiCV2Pred_perm, ngroups);
                % Permuted hold-out performance of the multi-class models
                % and its dichotomizers
                multicv2perf_perm               = nk_MultiClassAssessConfMatrix(MultiCV2ConfMatrix_perm, inp.labels, MultiCV2Pred_perm, MultiCV2Errs_perm, 'BAC');
                I.VCV2PERM_GLOBAL_MULTI(perms)  = multicv2perf_perm.BAC_Mean;
                I.VCV2PERM_GLOBAL_MULTI_ONEvsREST(:, perms) = multicv2perf_perm.BAC;
                % Evaluate permutation test criterion in multi-class model
                % and its dichotomizers
                crt_multi = feval(compfun, I.VCV2PERM_GLOBAL_MULTI(perms), I.VCV2MORIG_GLOBAL_MULTI.BAC_Mean);
                crt_onevsrest = feval(compfun, I.VCV2PERM_GLOBAL_MULTI_ONEvsREST(:,perms), I.VCV2MORIG_GLOBAL_MULTI.BAC');
                if ~crt_multi, fprintf('*'); else, fprintf('.'); end
                % Store test results
                I.VCV2MPERM_GLOBAL_MULTI(perms) = crt_multi;
                I.VCV2MPERM_GLOBAL_MULTI_ONEvsREST(:,perms) = crt_onevsrest;
            end
        end
    end
    
    for n=1:nM
        
        % Number of classifiers / predictors
        [ D, datatype, brainmaski, badcoordsi, labeli, labelopi ] = getD(FUSION.flag, inp, n);
        
        if iscell(inp.VIS), nVIS = inp.VIS{n}; else nVIS = inp.VIS; end
        
        % Probability of feature selection across all CV2 * CV1 partitions
        for h=1:nclass
            NumPredDiv = repmat(I.VCV2NMODEL(h),size(I.VCV2SEL{h,n},1),1);
            if ~decompfl 
                if PCV1SUMflag && h==1, I.PCV2 = zeros(D,nclass); end
                I.PCV2(:,h) = I.PCV2SUM{h, n}' ./ NumPredDiv; 
                if size(I.VCV2PEARSON{h, n},2)>1
                    % Compute mean and SD univariate association measures
                    % (currently Pearson and Spearman), and their uncorrected
                    % and FDR-corrected statistical significance
                    I.VCV2PEARSON_STD{h,n} = nm_nanstd(I.VCV2PEARSON{h, n},2);
                    I.VCV2PEARSON{h,n} = nm_nanmedian(I.VCV2PEARSON{h, n},2);
                    I.VCV2SPEARMAN_STD{h,n} = nm_nanstd(I.VCV2SPEARMAN{h, n},2);
                    I.VCV2SPEARMAN{h,n} = nm_nanmedian(I.VCV2SPEARMAN{h, n},2);
                    I.VCV2PEARSON_UNCORR_PVAL_STD{h,n} = nm_nanstd(I.VCV2PEARSON_UNCORR_PVAL{h, n},2);
                    I.VCV2PEARSON_UNCORR_PVAL{h,n} = nm_nanmedian(I.VCV2PEARSON_UNCORR_PVAL{h, n},2);
                    I.VCV2SPEARMAN_UNCORR_PVAL_STD{h,n} = nm_nanstd(I.VCV2SPEARMAN_UNCORR_PVAL{h, n},2);
                    I.VCV2SPEARMAN_UNCORR_PVAL{h,n} = nm_nanmedian(I.VCV2SPEARMAN_UNCORR_PVAL{h, n},2);
                    I.VCV2PEARSON_FDR_PVAL{h,n} = nm_nanmedian(I.VCV2PEARSON_FDR_PVAL{h, n},2);
                    I.VCV2SPEARMAN_FDR_PVAL{h,n} = nm_nanmedian(I.VCV2SPEARMAN_FDR_PVAL{h, n},2);
                    % Compute analytical p values using the method of
                    % Gaonkar et al. if linear SVM was the base learner
                    if linsvmfl
                        I.VCV2PVAL_ANALYTICAL{h,n} = nm_nanmedian(I.VCV2PVAL_ANALYTICAL{h, n},2);
                        I.VCV2PVAL_ANALYTICAL_FDR{h,n} = nm_nanmedian(I.VCV2PVAL_ANALYTICAL_FDR{h, n},2);
                    end
                end
            end
            % Compute empirical multivariate p values (corrected and
            % uncorrected as well as Z scores
            if any(permfl),
                I.VCV2ZSCORE{h,n}   = nm_nansum(I.VCV2ZSCORE{h,n},2) ./ I.VCV2SEL{h,n};
                % Uncorrected p values
                Pvals = nm_nansum(I.VCV2PERM{h,n},2)./ol ;
                Pvals(Pvals==0) = realmin;
                %Pvals    = 1-normcdf(I.VCV2ZSCORE{h,n}); 
                I.VCV2PERM{h,n} = -log10(Pvals);
                % FDR-Corrected p values
                Pvals = nm_nansum(I.VCV2PERM_FDR{h,n},2)./ol ;
                Pvals(Pvals==0) = realmin;
                I.VCV2PERM_FDR_PVAL{h,n} = -log10(Pvals);
            end
            
            % Mean relevance/weight vector across all CV2 * CV1 partitions
            I.VCV2{h,n} = nm_nansum(I.VCV2SUM{h,n},2) ./ I.VCV2SEL{h,n};
            
            % Compute standard error of relevance and weight vector 
            if size(I.VCV2SQ{h, n},2)>1
                I.VCV2SE{h,n}       = nm_nanmedian(sqrt(abs(I.VCV2SQ{h, n} - I.VCV2SUM2{h,n})./(I.VCV2SEL{h,n})),2);
                I.VCV2MEAN_CV1{h,n} = nm_nanmedian(I.VCV2MEAN{h,n},2);
                I.VCV2SE_CV1{h,n}   = nm_nanmedian(I.VCV2STD{h,n},2);
            else
                I.VCV2SE{h,n}  =     sqrt(abs(I.VCV2SQ{h, n} - I.VCV2SUM2{h,n})./(I.VCV2SEL{h,n}));
                I.VCV2MEAN_CV1{h,n} = I.VCV2MEAN{h,n};
                I.VCV2SE_CV1{h,n}   = I.VCV2STD{h,n};
            end
            
            % Compute CV-ratio 
            I.VCV2rat{h,n} = I.VCV2{h,n} ./ I.VCV2SE{h,n};
            
            % Compute Sign-based consistency, Z scores and P values
            [I.VCV2signconst_p{h,n}, I.VCV2signconst_pfdr{h,n}, I.VCV2signconst_z{h,n}, I.VCV2signconst{h,n}] = nk_SignBasedConsistencySignificance(I.VCV2VCV1{h,n});  
            
             % Compute GrandMean metrics
            I.VCV2rat_CV1{h,n}  = I.VCV2MEAN_CV1{h,n}./I.VCV2SE_CV1{h,n};
            I.VCV2MEANthreshSE_CV1{h,n} = zeros(size(I.VCV2MEAN_CV1{h,n}),'single');
            indMEANgrSE = abs(I.VCV2MEAN_CV1{h,n}) > I.VCV2SE_CV1{h,n};
            I.VCV2MEANthreshSE_CV1{h,n}(indMEANgrSE) = I.VCV2MEAN_CV1{h,n}(indMEANgrSE);

            % Compute voxel selection probability using 95% confidence interval
            % method
            I.VCV2PROB{h,n} = (nm_nansum(I.VCV2PROB{h,n},2)/ol).*sign(I.VCV2MEAN_CV1{h,n});
        end

        if numel(inp.tF)>1 && datatype
            fprintf('\n\nWriting out images for Modality #%g',i)
        end
        
        % Now we have to differentiate between imaging and non-imaging
        % analyses. In the former case we write out data to the disk
        if datatype ==1 || datatype==2
            
            currdir = pwd;
            cd(inp.rootdir);

            for h=1:nclass % Loop through binary classifiedars (predictors)
                % _______________________________________________________________________________
                % Generate filenames & save data:
                imgname = SAV.matname; 
                suff = ['_NumPred-' num2str(I.VCV2NMODEL(h))];
                varsuff = sprintf('_var%g',inp.tF(n));
                switch MODEFL
                    case 'regression'
                            basename ='PredictVol';
                            suff = [multlabelstr suff varsuff '_ID' id];
                    case 'classification'
                            basename = 'DiscrimVol';
                            suff = [multlabelstr '_cl' num2str(h) suff varsuff '_ID' id];
                end
                % _______________________________________________________________________________
                % Base images:
                vols = [I.VCV2{h,n} ...                 %  mean image
                        I.VCV2rat{h,n} ...              %  CV ratio image
                        I.VCV2SE{h,n} ...               %  standard error image
                        I.VCV2MEAN_CV1{h,n} ...         %  grand mean image (CV1 level)
                        I.VCV2SE_CV1{h,n} ...           %  grand mean standard error image (CV1 level)
                        I.VCV2MEANthreshSE_CV1{h,n} ... %  grand mean thresholded by (mean > standard error) image (CV1 level)
                        I.VCV2rat_CV1{h,n} ...          %  grand mean CV ratio image 
                        I.VCV2PROB{h,n} ...             %  grand mean CV ratio image 
                        I.VCV2signconst{h,n} ...        %  sign-based consistency image
                        I.VCV2signconst_z{h,n} ...      %  sign-based consistency Z score image
                        I.VCV2signconst_p{h,n} ...      %  sign-based consistency P value (uncorrected) image
                        I.VCV2signconst_pfdr{h,n} ...   %  sign-based consistency P value (FDR-corr) image
                        ];
                volnam = [        basename '_Mean_' imgname suff ];
                volnam = char(volnam, [basename '_CVratio_' imgname suff ]);
                volnam = char(volnam, [basename '_SE_' imgname suff] );
                volnam = char(volnam, [basename '_Mean-GrM_' imgname suff ]); 
                volnam = char(volnam, [basename '_SE-GrM_' imgname suff ]);
                volnam = char(volnam, [basename '_Mean-gr-SE-GrM_' imgname suff ]);
                volnam = char(volnam, [basename '_CVratio-GrM_' imgname suff ]);
                volnam = char(volnam, [basename '_Prob95CI-GrM_' imgname suff ]);
                volnam = char(volnam, [basename '_SignBased_' imgname suff ]);
                volnam = char(volnam, [basename '_SignBased_Z_' imgname suff ]);
                volnam = char(volnam, [basename '_SignBased_p_uncorr_' imgname suff ]);
                volnam = char(volnam, [basename '_SignBased_p_FDR_' imgname suff ]);
                % _______________________________________________________________________________
                % Univariate stats images:
                % Grand mean Spearman image 
                if exist('I.VCV2SPEARMAN','var') && ~isempty(I.VCV2SPEARMAN{h,n})
                    vols = [ vols I.VCV2SPEARMAN{h,n}];  
                    volnam = char(volnam, [basename '_Spearman-GrM_' imgname suff ]);
                end
                % Grand mean Spearman image 
                if exist('I.VCV2PEARSON','var') && ~isempty(I.VCV2PEARSON{h,n})
                    vols = [ vols I.VCV2PEARSON{h,n}];
                    volnam = char(volnam, [basename '_Pearson-GrM_' imgname suff ]);
                end
                % _______________________________________________________________________________
                % Permutation-based stats images
                if permfl 
                    vols = [vols I.VCV2ZSCORE{h,n} ...      % Save permutation-based Zscore image 
                                 I.VCV2PERM{h,n} ...         % Save permutation-based probability image(uncorrected)
                                 I.VCV2PERM_FDR_PVAL{h,n} ]; % Save permutation-based probability image (FDR-corrected)
                    volnam = char(volnam, [basename '_PermZ_' imgname suff ]);
                    volnam = char(volnam, [basename '_PermProb_' imgname suff ]);
                    volnam = char(volnam, [basename '_PermProbFDR_' imgname suff ]);
                end
                % _______________________________________________________________________________
                % Grand mean images:
                if isfield(nVIS,'mean')

                    mI.GCV2_t = nk_Threshold(mI.GCV2, inp.VIS{n}.mean); % Save grand mean image
                    sdI.GCV2_t = nk_Threshold(sdI.GCV2, inp.VIS{n}.se); % Save grand STD image
                    ind = ( mI.GCV2_t ~= 0 & sdI.GCV2_t > 0 );
                    normI.GCV2 = zeros(size(zI.GCV2)); 
                    normI.GCV2(ind) = zI.GCV2(ind);                     % Save grand Z score image
                    indCV2MEANgrCV2SE = abs(I.GCV2SUMh) > repmat(seI.GCV2,1,ol);
                    GrandProb95CI = sum(indCV2MEANgrCV2SE,2)./ol;       % Save grand probability image
                    vols = [vols mI.GCV2_t sdI.GCV2_t normI.GCV2 GrandProb95CI(:) ];    
                    volnam = char(volnam, [basename '_GrandMean_thresh_' imgname suff ]);
                    volnam = char(volnam, [basename '_GrandSTD_thresh_' imgname suff ]);
                    volnam = char(volnam, [basename '_GrandZ_thresh_' imgname suff ]);
                    volnam = char(volnam, [basename '_GrandProb95CI_' imgname suff ]);
                end

                clear mI.GCV2_t sdI.GCV2_t ind normI.GCV2 zI.GCV2 indCV2MEANgrCV2SE GrandProb95CI
                % _______________________________________________________________________________
                % Thresholded images:
                if isfield(nVIS,'thresh') && ~isempty(nVIS.thresh) && nVIS.thresh

                    I.VCV2_t    = nk_Threshold(I.VCV2{h,n}, inp.VIS{n}.mean);
                    I.VCV2SE_t  = nk_Threshold(I.VCV2SE{h,n}, inp.VIS{n}.se);
                    ind         = ( I.VCV2_t ~= 0 & I.VCV2SE_t > 0 );
                    I.VCV2rat_t = zeros(size(I.VCV2_t));
                    I.VCV2rat_t(ind) = I.VCV2_t(ind) ./ I.VCV2SE_t(ind);

                    vols = [vols I.VCV2_t I.VCV2rat_t I.VCV2SE_t ];

                    % Save thresholded mean image
                    threshstr = num2str(inp.VIS{n}.mean.val);
                    threshstr = regexprep(threshstr,' ','-');
                    volnam = char(volnam, [basename '_Mean_thresh-' threshstr '_' imgname suff ]);

                    % Save thresholded CV ratio image
                    volnam = char(volnam, [basename '_CVratio_thresh_' imgname suff ]);

                    % Save standard error image
                    threshstr = num2str(inp.VIS{n}.se.val);
                    threshstr = regexprep(threshstr,' ','-');
                    volnam = char(volnam, [basename '_SE_thresh-' threshstr '_' imgname suff ]);

                end 
                % _______________________________________________________________________________
                % Write-out:
                switch datatype
                    case 1 % SPM-based NIFTI write-out
                        nk_WriteVol(vols ,volnam, 2, brainmaski,[], labeli, labelopi);
                    case 2 % Surface-based write-out
                        s = MRIread(brainmaski);
                        for yy=1:size(vols,2)
                            filename = fullfile(pwd,[deblank(volnam(yy,:)) '.mgh']);
                            s.vol = vols(:,yy);
                            MRIwrite(s,filename)
                        end
                end
            end
            cd(currdir); 
        end
        
        %% Build output structure
        visdata{n}.params.dimvecx      = [1 D];
        visdata{n}.params.varind       = inp.tF(n);
        visdata{n}.params.visflag      = datatype;
        visdata{n}.params.brainmask    = brainmaski;
        visdata{n}.params.badcoords    = badcoordsi;
        visdata{n}.params.I.numCV2part = ll-1;
        visdata{n}.params.NumPred      = I.VCV2NMODEL(h);
        
        if ~isempty(featnames) && ~isempty(featnames{n});
            visdata{n}.params.features = featnames{n};
        else
            visdata{n}.params.features = cellstr(num2str((1:D)'));
        end
        visdata{n}.params.nfeats       = numel(visdata{n}.params.features);
        visdata{n}.MEAN                = I.VCV2(:,n);
        visdata{n}.SE                  = I.VCV2SE(:,n);
        visdata{n}.CVRatio             = I.VCV2rat(:,n);
        if isfield(I,'PCV2'), 
            visdata{n}.FeatProb        = {I.PCV2}; end
        visdata{n}.MEAN_CV2            = I.VCV2MEAN_CV1(:,n);
        visdata{n}.SE_CV2              = I.VCV2SE_CV1(:,n);
        visdata{n}.CVRatio_CV2         = I.VCV2rat_CV1(:,n);
        visdata{n}.Prob_CV2            = I.VCV2PROB(:,n);
        visdata{n}.SignBased_CV2       = I.VCV2signconst(:,n);
        visdata{n}.SignBased_CV2_p_uncorr = I.VCV2signconst_p(:,n);
        visdata{n}.SignBased_CV2_p_fdr = I.VCV2signconst_pfdr(:,n);
        visdata{n}.SignBased_CV2_z     = I.VCV2signconst_z(:,n);
        
        if ~decompfl(n)
            visdata{n}.Pearson_CV2              = I.VCV2PEARSON(:,n);
            visdata{n}.Spearman_CV2             = I.VCV2SPEARMAN(:,n);
            visdata{n}.Pearson_CV2_p_uncorr     = I.VCV2PEARSON_UNCORR_PVAL(:,n);
            visdata{n}.Spearman_CV2_p_uncorr    = I.VCV2SPEARMAN_UNCORR_PVAL(:,n);
            visdata{n}.Pearson_CV2_p_fdr        = I.VCV2PEARSON_FDR_PVAL(:,n);
            visdata{n}.Spearman_CV2_p_fdr       = I.VCV2SPEARMAN_FDR_PVAL(:,n);
            if isfield(I,'VCV2PEARSON_STD')
                visdata{n}.Pearson_CV2_STD          = I.VCV2PEARSON_STD(:,n);
                visdata{n}.Spearman_CV2_STD         = I.VCV2SPEARMAN_STD(:,n);
                visdata{n}.Pearson_CV2_p_uncorr_STD = I.VCV2PEARSON_UNCORR_PVAL_STD(:,n);
                visdata{n}.Spearman_CV2_p_uncorr_STD= I.VCV2SPEARMAN_UNCORR_PVAL_STD(:,n);
            end
            if linsvmfl
                visdata{n}.Analytical_p = I.VCV2PVAL_ANALYTICAL(:,n);
                visdata{n}.Analyitcal_p_fdr = I.VCV2PVAL_ANALYTICAL_FDR(:,n);
            end
        end
        if any(permfl)
            visdata{n}.PermProb_CV2             = I.VCV2PERM(:,n);
            visdata{n}.PermProb_CV2_FDR         = I.VCV2PERM_FDR(:,n); 
            visdata{n}.PermZ_CV2                = I.VCV2ZSCORE(:,n);
            visdata{n}.PermModel                = I.VCV2MODELP;
            visdata{n}.PermModel_std            = I.VCV2MODELP_STD;
            visdata{n}.PermModel_CV2            = I.VCV2MPERM_CV2;
            visdata{n}.PermModel_Eval_Global    = I.VCV2MPERM_GLOBAL;
            visdata{n}.PermModel_Crit_Global    = I.VCV2MPERM_EVALFUNC_GLOBAL;
            visdata{n}.ObsModel_Eval_Global     = I.VCV2MORIG_EVALFUNC_GLOBAL;
            if inp.multiflag
                visdata{n}.PermModel_Eval_Global_Multi_Bin  = I.VCV2MPERM_GLOBAL_MULTI_ONEvsREST;
                visdata{n}.PermModel_Crit_Global_Multi_Bin  = I.VCV2PERM_GLOBAL_MULTI_ONEvsREST;
                visdata{n}.ObsModel_Eval_Global_Multi_Bin   = I.VCV2MORIG_GLOBAL_MULTI.BAC;
                visdata{n}.PermModel_Eval_Global_Multi      = I.VCV2MPERM_GLOBAL_MULTI;
                visdata{n}.PermModel_Crit_Global_Multi      = I.VCV2PERM_GLOBAL_MULTI;
                visdata{n}.ObsModel_Eval_Global_Multi       = I.VCV2MORIG_GLOBAL_MULTI.BAC_Mean;
            end
        end
        if ~isempty(inp.extraL)
            visdata{n}.ExtraL = I.EXTRA_L;
        end
    end
end
