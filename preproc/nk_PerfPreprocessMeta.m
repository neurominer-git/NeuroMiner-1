function [tY, Pnt, paramfl, tYocv] = nk_PerfPreprocessMeta(inp, labels, paramfl)
% =========================================================================
% function [tY, Pnt, paramfl, tYocv] = nk_PerfPreprocessMeta(inp, labels, paramfl)
% =========================================================================
%
% OUTPUTS:
% -------
% tY        : generated decision-based data from lower-layer models
% Pnt       : preprocessing parameter structure containing information on
%              the stacking process performed
% paramfl   : the modified paramfl
% tYocv     : Independent test data
% 
% INPUTS:
% --------
% inp       : input data and some parameters
% labels    : the targets for prediction
% paramfl   : parameter structure describing the stacking process to be
%               done
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 01/2020

global MODEFL MULTI CV RAND VERBOSE PREPROC STACKING SVM

i = inp.f; j = inp.d; kbin = inp.nclass;

if ~exist('paramfl','var'), paramfl.use_exist = false; end
cv2flag = false; 
if isfield(PREPROC,'CV2flag') && (PREPROC.CV2flag - 1) == 1; cv2flag = true; end

% Binary/regression or multi-group processing mode
% if iscell(PREPROC)
%     BINMOD = PREPROC{1}.BINMOD;
% else
%     BINMOD = PREPROC.BINMOD;
% end
BINMOD = 1;

% Pointers
[iy,jy]     = size(CV.cvin{1,1}.TrainInd); 
TrInd       = CV.TrainInd{i,j}; 
TsInd       = CV.TestInd{i,j}; 

% Data containers
tY.Tr       = cell(iy,jy);
tY.CV       = cell(iy,jy);
tY.Ts       = cell(iy,jy);

% Labels & Indices
tY.TrL      = cell(iy,jy);
tY.CVL      = cell(iy,jy);
tY.TrInd    = cell(iy,jy);
tY.CVInd    = cell(iy,jy);
tY.VI       = cell(iy,jy);
% The order of statements here is critical: TsI should receive the original
% TsInd before it is modified below!
%
% We have to ignore BINMOD for multi-class stacking operating on
% dichotomization-based input channels
inp.multiflag=0;
if ~isempty(MULTI) && MULTI.flag, tY.mTsL = labels(TsInd,:); end

for u=1:kbin
    switch MODEFL
        case 'regression'
            TsL         = labels(TsInd,:);
            TsInd       = true(size(TsL,1),1);
        case 'classification'
            if RAND.Decompose == 9
                TsL     = labels(TsInd,:);
                TsInd   = true(size(TsL,1),1);
            else
                TsL     = CV.classnew{i,j}{u}.label;
                TsInd   = CV.classnew{i,j}{u}.ind;
            end
    end
    tY.TsL{u} = TsL;
    tY.TsInd{u} = TsInd;
end

ll = 0;
ukbin = kbin;   
if VERBOSE; fprintf('\nProcessing Mode: only binary / regression preprocessing allowed for stacking'); end
TsI = cell(kbin,1);
for u=1:ukbin, TsI{u} = TsInd; end

% Create preprocessing parameter structure if it has not been provided to
% the function
if ~exist('paramfl','var') || ~paramfl.found || ~paramfl.use_exist
    Pnts = struct('data_ind', [], ...
                  'train_ind', [], ...
                  'nP', [], ...
                  'nA', [], ...
                  'TrainedParam', []);
    Pnt = repmat(Pnts,iy,jy,ukbin);
else
    Pnt = paramfl.Param;
    paramfl.Param = [];
end

% Number of input channels / analyses
nA = numel(inp.analyses); cnt=1; flfnd=false; 

% Model count per modality
tY.nM_cnt = zeros(1,nA);

% Input data container gathering CVdatamats across analyses
GDD = struct('RootPath',[],'CVfilename',[],'GD',[]);

% Check whether the function is running in OOCV mode and determine whether
% the OOCV module has been for each and every input channel of the stacker
if isfield(inp,'oocvind') && ~isempty(inp.oocvind)
    oocvflag = true;
    oocvind = inp.oocvind;
    oocvconducted = false(1,nA);
    for am = 1:nA
        if isfield(inp.analyses{am},'OOCV')
            oocvconducted(am) = true;
        end
    end
    smiss = sum(oocvconducted==0); missf = find(oocvconducted==0);
    if smiss>0
       fprintf('\n');cprintf('red*','=================================================================================================== ')
       fprintf('\n'); cprintf('red*','%g input-layer model(s) have not been applied to the test data! The stacked model needs all of them. ', smiss);
       fprintf('\n'); cprintf('red*','Following models need to be applied to the independent data to generate features for the stacker: ')
       for amx = 1:smiss
           fprintf('\n'); cprintf('red*','==> Analysis ID: %s ', inp.analyses{missf(amx)}.id);
       end
       error('Aborting stacker prepocessing module due to missing input features')
    else
        % Container for OOCV data
        OOCVDD = struct('RootPath',[],'FileNames',[], 'GD', [], 'CNT', []);
    end
else
    oocvflag = false;
end

if oocvflag, tYocv.Ts = cell(iy,jy); else, tYocv = []; end

% Get analysis IDs
analid = cell(1, nA); for am = 1:nA, analid{am} = inp.analyses{am}.id; end
analid=char(analid);

% Here we collect the filepaths to the CVdatamats [and OOCVdatamats, if OOCV mode is used]
for am = 1:nA
    
    nM = numel(inp.analyses{am}.GDdims);
    % Loop through modalities (for fusion settings)
    for jm=1:nM
        if nM>1, bgstr = 'bagged '; else, bgstr=''; end
         if ~exist(inp.analyses{am}.GDdims{jm}.RootPath,'dir')
             cprintf('red*','\nRoot path %s of analysis %g not found.', inp.analyses{am}.GDdims{jm}.RootPath, am);
             GDD(cnt).RootPath = cellstr(spm_select(nM,'dir',sprintf('Select the root path(s) [n=%g] of %sanalysis %g [%s]:', nM, bgstr, am, inp.analyses{am}.desc)));
         else
            GDD(cnt).RootPath = inp.analyses{am}.GDdims{jm}.RootPath;
         end
         GDD(cnt).CVfilename = inp.analyses{am}.GDdims{jm}.GDfilenames{inp.f,inp.d};
         if iscell(GDD(cnt).RootPath)
             for xu=1:numel(RootPath)
                CVpath = fullfile(GDD(cnt).RootPath{xu}, GDD(cnt).CVfilename);
                if exist(CVpath,'file'), flfnd = true; break; end
             end
         else
             CVpath = fullfile(GDD(cnt).RootPath,GDD(cnt).CVfilename);
             if exist(CVpath,'file'), flfnd = true; end
         end
         if ~flfnd || exist(CVpath,'dir'),
             error('CVdatamat is missing for:\tAnalysis %g (%s),\n\tModality %g,\n\tCV2 [%g, %g].\nMake sure you have computed it before running the stacked analysis.', am, inp.analyses{am}.id, jm, inp.f, inp.d);
         end
         cprintf('blue*','\n%s: ', analid(am,:)); [~,nam]= fileparts(CVpath); fprintf('Loading %s', nam); load(CVpath,'GD'); 
         GDD(cnt).GD = GD; clear GD;
         
         if oocvflag
            OOCVDD(cnt).RootPath    = inp.analyses{am}.OOCV{oocvind}.RootPath{jm};
            OOCVDD(cnt).FileNames   = inp.analyses{am}.OOCV{oocvind}.FileNames{jm}{inp.f,inp.d};
            OOCVpath                = fullfile(OOCVDD(cnt).RootPath,[OOCVDD(cnt).FileNames, '.mat']);
            if ~exist(OOCVpath,'file')
                error('OOCVdatamat is missing for Analysis %g, Modality %g, CV2 [%g, %g].\nMake sure you have computed it before running the stacked analysis.',  am, jm, inp.f, inp.d);
            end
            load(OOCVpath,'binOOCVDh','cntOOCVDh'); 
            OOCVDD(cnt).GD          = binOOCVDh; 
            OOCVDD(cnt).CNT         = cntOOCVDh; 
         end
         cnt=cnt+1;
    end
end
fprintf('\nDone loading CVdatamats.')

% train parameters eventually only for a subset of partitions
if isfield(inp,'CV1p')
    sta_iy = inp.CV1p(1); stp_iy = inp.CV1p(2);
    sta_jy = inp.CV1f(1); stp_jy = inp.CV1f(2);
else
    sta_iy = 1; stp_iy = iy;
    sta_jy = 1; stp_jy = jy;
end

for k=sta_iy:stp_iy % Inner permutation loop

    for l=sta_jy:stp_jy % Inner CV fold loop

        ll = ll+1;
        tElapsed = tic;
        fprintf('\nWorking on CV1 [ %2g, %2g ]: Prepare data', k, l);
        
        if isempty(tY.Tr{k,l})
            tY.Tr{k,l} = cell(1,ukbin); tY.CV{k,l} = cell(1,ukbin); tY.Ts{k,l} = cell(1,ukbin); tY.VI{k,l} = cell(1,ukbin);
        end
       
        for u=1:ukbin % Binary comparison loop depending on PREPROC.BINMOD 

            switch MODEFL
                case 'classification'
                    if BINMOD == 1 && ~MULTI.flag 
                            I = CV.class{i,j}{u};
                    else
                            I = CV.cvin{i,j};
                    end
                case 'regression'
                    I = CV.cvin{i,j};
            end
            TrI = TrInd(I.TrainInd{k,l});
            CVI = TrInd(I.TestInd{k,l});
            TrL = labels(TrI,:);
            CVL = labels(CVI,:);

            if size(TrI,2)>2, TrI = TrI'; CVI = CVI'; TsI{u} = TsI{u}'; end
            trd = []; cvd=[]; tsd = []; ocv = []; 
            cnt=1; chn = []; chn_cnt = 0;
            
            for am = 1:nA

                nM = numel(inp.analyses{am}.GDdims);
                
                for jm = 1:nM
                    
                    [~,Pspos, nP] = nk_GetModelParams2(inp.analyses{am}.GDdims{jm}, inp.multiflag, inp.ll, u);
                    CVfilename = inp.analyses{am}.GDdims{jm}.GDfilenames{inp.f,inp.d}; 
                    tY.nM_cnt(am) = nP;
                    if VERBOSE, 
                        fprintf('\nCV [ %g, %g ]: Concatenating predictions from modality #%g in analysis #%g (%s): %s',k, l,jm, am, inp.analyses{am}.id, CVfilename); 
                        if oocvflag
                            OOCVfilename = inp.analyses{am}.OOCV{oocvind}.FileNames{jm}{inp.f,inp.d};
                            fprintf('\nOOCV [ %g, %g ]: Concatenating predictions from modality #%g in analysis #%g (%s): %s',k, l,jm, am, inp.analyses{am}.id, OOCVfilename); 
                        end
                    else
                        fprintf('.');
                    end
                    % Load CVdatamat
                    chn_cnt = chn_cnt + 1;
                    i_trd = []; i_cvd = []; i_tsd = []; i_ocv = [];
                    
                    for n=1:numel(Pspos)

                        % for stacked generalization we need to find the
                        % out-of-training prediction across the
                        % crossvalidation cycle
                        mdfl=false;
                        switch STACKING.mode
                            case 1
                                if size(GDD(cnt).GD.DV{Pspos(n)}{k,l,u},2)>1
                                   mdfl = true;
                                   i_trd = [i_trd nm_nanmedian(double(GDD(cnt).GD.DT{Pspos(n)}{k,l,u}))];
                                else
                                   i_trd = [i_trd double(GDD(cnt).GD.DT{Pspos(n)}{k,l,u})]; 
                                end
                            case 2
                                oot_ind  = [];
                                oot_data = [];
                                % First build an index array and data array
                                % of out of fold predictions
                                for lll = 1:jy
                                   oot_ind = [oot_ind; I.TestInd{k,lll}];
                                   if size(GDD(cnt).GD.DV{Pspos(n)}{k,lll,u},2)>1
                                       mdfl = true;
                                       oot_data = [oot_data; nm_nanmedian(double(GDD(cnt).GD.DV{Pspos(n)}{k,lll,u}),2)];
                                   else
                                       oot_data = [oot_data; double(GDD(cnt).GD.DV{Pspos(n)}{k,lll,u})];
                                   end
                                end
                                [~,oot_order] = ismember(I.TrainInd{k,l},oot_ind);
                                i0 = oot_order == 0;
                                if any(i0), oot_order(i0) = I.TrainInd{k,l}(i0); end
                                i_trd = [i_trd oot_data(oot_order,:)];
                        end
                        if size(GDD(cnt).GD.DV{Pspos(n)}{k,l,u},2)>1
                            i_cvd = [i_cvd nm_nanmedian(double(GDD(cnt).GD.DV{Pspos(n)}{k,l,u}),2)];
                            i_tsd = [i_tsd nm_nanmedian(double(GDD(cnt).GD.DS{Pspos(n)}{k,l,u}),2)];
                            if oocvflag, 
                                Is = OOCVDD(cnt).CNT{u}{n}(k,l,1);
                                Ie = OOCVDD(cnt).CNT{u}{n}(k,l,2);
                                i_ocv = [ i_ocv nm_nanmedian(double(OOCVDD(cnt).GD{u}(:,Is:Ie)),2)]; 
                            end
                        else
                            i_cvd = [i_cvd double(GDD(cnt).GD.DV{Pspos(n)}{k,l,u})];
                            i_tsd = [i_tsd double(GDD(cnt).GD.DS{Pspos(n)}{k,l,u})];
                            if oocvflag, 
                                Is = OOCVDD(cnt).CNT{u}{n}(k,l,1);
                                Ie = OOCVDD(cnt).CNT{u}{n}(k,l,2);
                                i_ocv = [ i_ocv double(OOCVDD(cnt).GD{u}(:,Is:Ie))]; 
                            end
                        end   
                    end
                    if strcmp(SVM.prog,'SEQOPT') && numel(Pspos) > 1
                        numNode = 1;
                        i_trd = nm_nanmedian(i_trd,2);
                        i_cvd = nm_nanmedian(i_cvd,2);
                        i_tsd = nm_nanmedian(i_tsd,2);
                        if oocvflag, i_ocv = nm_nanmedian(i_ocv,2); end
                    else
                        numNode = size(i_trd,2);
                    end
                    trd = [trd i_trd];
                    cvd = [cvd i_cvd];
                    tsd = [tsd i_tsd];
                    if oocvflag, ocv = [ocv i_ocv] ; end
                    chn = [chn repmat(chn_cnt,1,numNode)];
                    cnt = cnt + 1;
                end
            end
            clear cnt    
            [InputParam.Tr{1},~, SrcParam.iTrX]     = nk_ManageNanCases(trd, TrL);
            [InputParam.Ts{1}, TrL, SrcParam.iTr]   = nk_ManageNanCases(trd, TrL);
            [InputParam.Ts{2}, CVL, SrcParam.iCV]   = nk_ManageNanCases(cvd, CVL);
            [InputParam.Ts{3}, ~, SrcParam.iTs]     = nk_ManageNanCases(tsd);
                    
            if oocvflag , InputParam.Ts{4} = ocv; end

            %% Generate SrcParam structure
            SrcParam.TrX                = TrI;
            SrcParam.TrI                = TrI;
            SrcParam.CVI                = CVI;
            SrcParam.TsI                = TsI{u};
            SrcParam.u                  = u;
            SrcParam.binmult            = 1;
            SrcParam.CV1perm            = k;
            SrcParam.CV1fold            = l;
            SrcParam.covars             = inp.covars;
            SrcParam.covars_oocv        = inp.covars_oocv;
            %SrcParam.orig_labels        = inp.labels;
            
            %% Run ADASYN if needed
            if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
                if VERBOSE, fprintf('\nUsing ADASYN to generate synthetic training data for partition CV2 [%g, %g], CV [%g, %g]', i, j, k, l); else; fprintf('\t...ADASYN'); end
                Covs = [];
                % Do we have covars? if so, they have to be integrated
                % into the creation of synthetic data.
                if ~isempty(SrcParam.covars)
                    Covs = SrcParam.covars( SrcParam.TrX(~SrcParam.iTrX),:);
                    % Check if covariance matrix contains NaNs and
                    % impute them.
                    if any(isnan(Covs(:)))
                        IN = struct('method','seuclidean','k',7,'X',Covs);
                        Covs = nk_PerfImputeObj(Covs, IN);
                    end
                end
                [ vTrSyn, LabelSyn, CovarsSyn ] = nk_PerfADASYN( InputParam.Tr, TrL, SVM.ADASYN, Covs, true);
                InputParam.TrSyn = vTrSyn;
                SrcParam.TrainLabelSyn = LabelSyn;
                if ~isempty(SrcParam.covars), SrcParam.covarsSyn = CovarsSyn; end
            end

            switch MODEFL
                case 'classification'
                    if RAND.Decompose ~=9
                        SrcParam.BinaryTrainLabel   = CV.class{i,j}{u}.TrainLabel{k,l};
                        SrcParam.BinaryCVLabel      = CV.class{i,j}{u}.TestLabel{k,l};
                    end
                    SrcParam.MultiTrainLabel    = TrL;
                    SrcParam.MultiCVLabel       = CVL;
                case 'regression'
                    SrcParam.TrainLabel         = TrL;
                    SrcParam.CVLabel            = CVL;
            end

            %% Generate and execute for given CV1 partition preprocessing sequence
            [InputParam, oTrainedParam, SrcParam] = nk_GenPreprocSequence(InputParam, PREPROC, SrcParam, Pnt(k, l, u).TrainedParam);
            
             % Check whether we have imputed labels
            if isfield(SrcParam,'TrL_imputed'), 
                TrL = SrcParam.TrL_imputed; 
                [~,TrL] = nk_ManageNanCases(InputParam.Ts{1}, TrL, SrcParam.iTr); 
                tTrL = labels(TrI,:);
                TrL(~isnan(tTrL)) = tTrL(~isnan(tTrL));
            else
                % Overwrite labels adjusted to NaN cases
                TrL = labels(TrI(~SrcParam.iTr),:); CVL = labels(CVI(~SrcParam.iCV),:);
            end
            
            %% Write preprocessed data to mapY structure
            if isfield(paramfl,'PREPROC') && isfield(paramfl,'PXfull') && ~isempty(paramfl.PXopt{u})
                % Here an optimized parameter space exists that has
                % been used to limit the computation to the unique
                % parameter cominations. We create and store the pointers
                % and used them later on to retrieve the right preproc
                % version of the data and the respective preproc params
                [ Pnt(k,l,u).data_ind, ...
                  Pnt(k,l,u).train_ind, ...
                  Pnt(k,l,u).nP, ...
                  Pnt(k,l,u).nA] = nk_ParamReplicator(paramfl.P{u}, paramfl.PXopt{u}, paramfl.PREPROC, oTrainedParam);
            end
            trd = InputParam.Ts(:,1); 
            cvd = InputParam.Ts(:,2); 
            tsd = InputParam.Ts(:,3);
            if oocvflag, ocv = InputParam.Ts(:,4); end
            
            if  any(cellfun(@isempty,trd)) ||  any(cellfun(@isempty,cvd)) ||  any(cellfun(@isempty,tsd)),
                error('Empty training/test/validation data matrices return by preprocessing pipeline. Check your settings')
            end
            
            switch BINMOD

                case 0 % Multi-group mode both in FBINMOD and PREPROC.BINMOD
                    if ukbin > 1
                        [tY.Tr{k,l}{u},TrL] = nk_ManageNanCases(trd, TrL, SrcParam.iTr);
                        [tY.CV{k,l}{u},CVL] = nk_ManageNanCases(cvd, CVL, SrcParam.iCV);
                        tY.Ts{k,l}{u} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                        if oocvflag, tYocv.Ts{k,l}{u} = ocv; end
                    else
                        [tY.Tr{k,l},TrL] = nk_ManageNanCases(trd, TrL, SrcParam.iTr);
                        [tY.CV{k,l},CVL] = nk_ManageNanCases(cvd, CVL, SrcParam.iCV);
                        tY.Ts{k,l} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                        if oocvflag, tYocv.Ts{k,l} = ocv; end
                    end
            
                    % Write dichotomization labels to CV1 partition
                    for zu=1:kbin % Binary loop depending on the # of binary comparisons
                        % Generate logical indices
                        if isfield(CV,'class') && length(CV.class{i,j}{zu}.groups) == 2
                            % One-vs-One
                            indtr = ( TrL == CV.class{i,j}{zu}.groups(1) | TrL == CV.class{i,j}{zu}.groups(2) ) | isnan(TrL);   
                            indcv = ( CVL == CV.class{i,j}{zu}.groups(1) | CVL == CV.class{i,j}{zu}.groups(2) ) | isnan(CVL);
                        else
                            % One-vs-All
                            indtr = TrL~=0;   
                            indcv = CVL~=0;
                        end
                        % Write indices
                        tY.TrInd{k,l}{zu} = indtr;
                        tY.CVInd{k,l}{zu} = indcv;
                        % Write labels to CV1 partition
                        if isfield(CV,'class')
                            tY.TrL{k,l}{zu} = CV.class{i,j}{zu}.TrainLabel{k,l};
                            tY.CVL{k,l}{zu} = CV.class{i,j}{zu}.TestLabel{k,l};	
                        else
                            tY.TrL{k,l}{zu} = labels(CV.TrainInd{i,j}(CV.cvin{i,j}.TrainInd{k,l}));
                            tY.CVL{k,l}{zu} = labels(CV.TrainInd{i,j}(CV.cvin{i,j}.TestInd{k,l}));
                        end
                        tY.VI{k,l}{zu} = chn;
                    end
                case 1
                    % Write data to CV1 partition
                    [tY.Tr{k,l}{u},TrL] = nk_ManageNanCases(trd, TrL, SrcParam.iTr); 
                    [tY.CV{k,l}{u},CVL] = nk_ManageNanCases(cvd, CVL, SrcParam.iCV); 
                    [tY.Ts{k,l}{u}] = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                    if oocvflag, tYocv.Ts{k,l}{u} = ocv; end
                    
                    if ~strcmp(MODEFL,'regression') && length(CV.class{i,j}{u}.groups) == 2
                        indtr = (TrL == CV.class{i,j}{u}.groups(1) | TrL == CV.class{i,j}{u}.groups(2)) | isnan(TrL);   
                        indcv = (CVL == CV.class{i,j}{u}.groups(1) | CVL == CV.class{i,j}{u}.groups(2)) | isnan(CVL);
                    else
                        indtr = true(size(TrI,1),1);   
                        indcv = true(size(CVI,1),1);
                    end
                    % Write indices
                    tY.TrInd{k,l}{u} = indtr;
                    tY.CVInd{k,l}{u} = indcv;

                    switch MODEFL
                        case 'regression' 
                            tY.TrL{k,l}{u} = TrL;
                            tY.CVL{k,l}{u} = CVL;
                        case 'classification'
                            if RAND.Decompose ~=9
                                % Write labels to CV1 partition
                                tY.TrL{k,l}{u} = CV.class{i,j}{u}.TrainLabel{k,l};	
                                tY.CVL{k,l}{u} = CV.class{i,j}{u}.TestLabel{k,l};	
                            else
                                tY.TrL{k,l}{u} = TrL;
                                tY.CVL{k,l}{u} = labels(CV.TrainInd{i,j}(CV.cvin{i,j}.TestInd{k,l}),:);
                            end
                    end
                    for jj=1:size(trd,1), tY.VI{k,l}{u}{jj} = chn'; end
            end

            if ~isempty(MULTI) && MULTI.flag, tY.mTrL{k,l} = TrL; tY.mCVL{k,l} = CVL; end    
            
            if paramfl.write || cv2flag, Pnt(k,l,u).TrainedParam = oTrainedParam; end
            clear InputParam TrainParam SrcParam
        end
        fprintf('\tCompleted in %1.2fs.',toc(tElapsed)); 
    end
end

%% Optionally, perform information integration across CV1 folds
% ... to be completed
