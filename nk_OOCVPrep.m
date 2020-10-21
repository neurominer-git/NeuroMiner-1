 function [ act, inp ] = nk_OOCVPrep(act, inp, parentstr)
% =========================================================================
% FORMAT [act, inp] = nk_OOCVprep(act, inp, parentstr)
% =========================================================================
% Runtime model of independent test data prediction
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, last modified 01/2020

global CV NM

% Detect completed analyses
complvec = []; for z=1:numel(NM.analysis), if NM.analysis{z}.status, complvec = [ complvec z ]; end; end

% Initialize runtime parameters 
if ~exist('inp','var') || isempty(inp)
    inp = struct( 'analind', complvec(1), ...
                    'oocvind', 1, ...
                    'OO', NM.OOCV{1}, ...
                    'lfl', 1, ...
                    'ovrwrt', 2, ...
                    'saveparam', 2, ...
                    'saveCV1', 2, ...
                    'loadparam', 2, ...
                    'batchflag', 0);
end
na_str = '?'; inp.datatype = 'OOCVdatamat'; 
OverWriteStr = []; GridSelectStr = []; LoadModelsStr = []; LoadParamsStr = []; LoadStr = []; SaveStr = []; SaveCV1Str = [];
OverWriteAct = []; GridSelectAct = []; LoadModelsAct = []; LoadParamsAct = []; LoadAct = []; SaveAct = []; SaveCV1Act = [];
DATASCRAM = false; if isfield(NM.defs,'data_scrambled') && ~isempty(NM.defs.data_scrambled), DATASCRAM = NM.defs.data_scrambled;end
    
%% Configure menu
% Select analysis
if numel(inp.analind)<2,
    AnalSelStr = sprintf('Analysis %g', inp.analind);   
else
    AnalSelStr = sprintf('%g Analyses: %s',numel(inp.analind), strjoin(cellstr(num2str(inp.analind'))',', '));
end
AnalSelectStr = sprintf('Choose analysis to work on [ %s ]|', AnalSelStr);                                                  AnalSelectAct = 1;   
analysis      = NM.analysis{inp.analind(1)}; 
% Select independent test data container
if isfield(inp,'oocvind'), OOCVSelStr = sprintf('New data #%g: %s', inp.oocvind, inp.OO.desc); else, OOCVSelStr = na_str; end
OOCVSelectStr = sprintf('Choose independent data to work on [ %s ]|', OOCVSelStr);                                          OOCVSelectAct = 2;   
if DATASCRAM, inp.loadparam = 1; inp.saveparam = 2; end
if ~isempty(analysis)
    
    % Initialize global parameters for the selected analysis
    nk_SetupGlobVars2(analysis.params, 'setup_main', 0); 
    
    % Compute from scratch or use pre-computed datamats ?
    LFL_opts        = {'Compute from scratch',sprintf('Use precomputed %s',inp.datatype)};                                      
    ModeStr         = sprintf('Operation mode of independent test module [ %s ]|',LFL_opts{inp.lfl});                       ModeAct = 3;
    
    if inp.lfl == 1,
        % from scratch
        OVRWRT_opts     = {'Overwrite existing','Do not overwrite'};       
        OverWriteStr = sprintf('Overwrite existing %s files [ %s ]|', inp.datatype, OVRWRT_opts{inp.ovrwrt}) ;              OverWriteAct = 4; 
    else
        % precomputed
        nOOCVFiles = na_str; 
        if isfield(inp,'oocvmat') && ~isempty(inp.oocvmat), 
            selGrid = ~cellfun(@isempty,inp.oocvmat); inp.GridAct = selGrid;
            nOOCVFiles = sprintf('%g selected', sum(selGrid(:))); 
        end     
        OverWriteStr = sprintf('Specify %s files [ %s ]|', inp.datatype, nOOCVFiles);                                       OverWriteAct = 4; 
    end
    
    % Retrieve CV2 partitions to operate on
    if ~isfield(inp,'GridAct'), inp.GridAct = analysis.GDdims{1}.GridAct; end;                                              
    GridSelectStr = sprintf('Select CV2 partitions to operate on [ %g selected ]|',  sum(inp.GridAct(:)));                  GridSelectAct = 5;
    
    % Configure loading of pre-existing parameters and models
    if inp.saveparam == 2 && inp.lfl == 1
        LOAD_opts        = {'yes', 'no'}; 
        if ~DATASCRAM
            LoadStr = sprintf('Use saved pre-processing params and models [ %s ]|', LOAD_opts{inp.loadparam});              LoadAct = 7;
        end
        if inp.loadparam == 1
            if isfield(inp,'optpreprocmat'), 
                selGridPreproc = ~cellfun(@isempty,inp.optpreprocmat);
                nParamFiles = sprintf('%g files selected', sum(selGridPreproc(:))); 
            else, 
                nParamFiles = na_str; 
            end
            LoadParamsStr = sprintf('Select preprocessing parameter files [ %s ]|' ,nParamFiles);                           LoadParamsAct = 8;
            if isfield(inp,'optmodelmat'), 
                selGridModel = ~cellfun(@isempty,inp.optmodelmat);
                nModelFiles = sprintf('%g files selected', sum(selGridModel(:))); 
            else, 
                nModelFiles = na_str; 
            end
            LoadModelsStr = sprintf('Select model files [ %s ]|',nModelFiles);                                              LoadModelsAct = 9;
        end
    end
    
    % If loading of pre-existing models and params is not chosen, allow to
    % save the computed params and models to disk
    if inp.loadparam == 2 && inp.lfl == 1
        SAVE_opts       = {'yes', 'no'};   
        SaveStr = sprintf('Save pre-processing params and models to disk [ %s ]|', SAVE_opts{inp.saveparam});               SaveAct = 6;
        if inp.saveparam == 1
            SaveCV1Str = sprintf('Save pre-processing params at CV1 level [ %s ]|', SAVE_opts{inp.saveCV1});                SaveCV1Act = 12;
        end
    end
end
 
 %% Build interactive menu
menustr = [ AnalSelectStr ...
            OOCVSelectStr ...
            ModeStr ...
            OverWriteStr ...
            GridSelectStr ...
            SaveStr ...
            SaveCV1Str ...
            LoadStr ...
            LoadParamsStr ... 
            LoadModelsStr ];

menuact = [ AnalSelectAct ...
            OOCVSelectAct ...
            ModeAct ...
            OverWriteAct ...
            GridSelectAct ...
            SaveAct ...
            SaveCV1Act ...
            LoadAct ...
            LoadParamsAct ...
            LoadModelsAct ];       

disallow = false;

%% Check whether all parameters are available
if ~sum(inp.GridAct(:)) || isempty(inp.analind), disallow = true; end

if inp.loadparam == 1
    if ~isfield(inp,'optpreprocmat') || isempty(inp.optpreprocmat), disallow = true; end
    if ~isfield(inp,'optmodelmat') || isempty(inp.optmodelmat), disallow = true; end
end

if ~disallow, menustr = [menustr '|PROCEED >>>']; menuact = [menuact 10]; end 

%% Display menu and act on user selections
nk_PrintLogo
mestr = 'Independent test module run-time configuration'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr);
if ~inp.batchflag && act<13, act = nk_input(mestr, 0, 'mq', menustr, menuact); end

switch act
    
    case 0
        return
    % Select analysis
    case 1 
        showmodalvec = []; analind = inp.analind; 
        if length(NM.analysis)>1, t_act = 1; brief = 1;
            while t_act>0, 
                [t_act, analind, ~, showmodalvec , brief] = nk_SelectAnalysis(NM, 0, navistr, analind, [], 1, showmodalvec, brief); 
            end; 
            if ~isempty(analind), inp.analind = complvec(analind) ; end
        end
        inp.GridAct = NM.analysis{inp.analind(1)}.GDdims{1}.GridAct;
    % Select OOCV data
    case 2    
        [ NM, OO, oocvind ] = nk_SelectOOCVdata(NM, 1, 0);  
        if ~isempty(oocvind), inp.OO = OO; inp.oocvind = oocvind; end 
    case 3
         lfl = nk_input('Define run-time mode of independent test module',0,'mq',strjoin(LFL_opts, '|'),[1,2],inp.lfl);
         if lfl, inp.lfl = lfl; end 
    % Overwrite?
    case 4
        switch inp.lfl
            case 1
                if inp.ovrwrt == 1, inp.ovrwrt = 2; elseif inp.ovrwrt  == 2, inp.ovrwrt = 1; end
            case 2
                tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind);
                inp.oocvmat = nk_GenDataMaster(NM.id, 'OOCVdatamat', CV, [], tdir);
        end
    case 5
        [operms,ofolds] = size(CV.TrainInd);
        tact = 1; while tact > 0 && tact < 10, [ tact, inp.GridAct ] = nk_CVGridSelector(operms, ofolds, inp.GridAct, 0); end
    case 6
        if inp.saveparam == 1, inp.saveparam = 2; elseif inp.saveparam == 2,  inp.saveparam = 1; end
    case 7
        if inp.loadparam == 1, inp.loadparam = 2; elseif inp.loadparam == 2,  inp.loadparam = 1; end
    case 8
        tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind);
        optpreprocmat = nk_GenDataMaster(NM.id, 'OptPreprocParam', CV, [], tdir);
        if ~isempty(optpreprocmat), inp.optpreprocmat = optpreprocmat; end
    case 9
        tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind);
        optmodelmat = nk_GenDataMaster(NM.id, 'OptModel', CV, [], tdir);
        if ~isempty(optmodelmat), inp.optmodelmat = optmodelmat; end
    case {10,11}
         inp.oocvname = sprintf('OOCV_%g',inp.oocvind);
         nA = 1; if numel(inp.analind)>1, nA = numel(inp.analind); end
         for i=1:nA
            nk_SetupGlobVars2(NM.analysis{inp.analind(i)}.params, 'setup_main', 0); 
            NM.runtime.curanal = inp.analind(i);
            inp.analysis_id = NM.analysis{inp.analind(i)}.id;
            inp.saveoptdir = [ NM.analysis{inp.analind(i)}.rootdir filesep 'opt' ];
            NM.analysis{inp.analind(i)}.OOCV{inp.oocvind} = OOCVPrep(NM, inp, NM.analysis{inp.analind(i)});
            nk_SetupGlobVars2(NM.analysis{inp.analind(i)}.params, 'clear', 0); 
         end
    case 12
        if inp.saveCV1 == 1, inp.saveCV1 = 2; elseif inp.saveCV1 == 2,  inp.saveCV1 = 1; end
end

function tdir = create_defpath(analysis, oocvind)
 
rootdir = analysis.rootdir;
algostr = getAlgoStr(analysis);
if isfield(analysis,'OOCV') && numel(analysis.OOCV) >= oocvind && isfield(analysis.OOCV{oocvind},'RootPath')
    if iscell(analysis.OOCV{oocvind}.RootPath)
        tdir = analysis.OOCV{oocvind}.RootPath{1};
    else
        tdir = analysis.OOCV{oocvind}.RootPath;
    end
else
    oocvdir = sprintf('OOCV_%g', oocvind);
    tdir = fullfile(rootdir, algostr, oocvdir);
end
%
% =========================================================================
function OOCVres = OOCVPrep(dat, inp1, analysis)
global SAV MODEFL CV OOCV FUSION MULTILABEL

if inp1.saveparam   == 2, inp1.saveparam    = 0; end
if inp1.loadparam   == 2, inp1.loadparam    = 0; end
if inp1.ovrwrt      == 2, inp1.ovrwrt       = 0; end
if inp1.lfl         == 1, inp1.analmode     = 0; else inp1.analmode = 1; end

F = 1; nF = 1;
if ~isempty(FUSION)        
    F = analysis.params.TrainParam.FUSION.M;
    nF = numel(F); if FUSION.flag < 3, nF = 1; end
    inp1.nF = nF;
end

if strcmp(MODEFL,'classification')
    inp1.nclass = length(CV.class{1,1});
else
    inp1.nclass = 1;
end

if isfield(inp1.OO,'label') && ~isempty(inp1.OO.label), 
    inp1.LabelCV     = dat.label; 
    inp1.labelOOCV   = inp1.OO.label; 
end
inp1.cases_oocv      = inp1.OO.cases;
inp1.nOOCVsubj       = numel(inp1.OO.cases);
inp1.id              = dat.id;
stranalysis          = SAV.matname;
inp1.ngroups         = numel(unique(dat.label(~isnan(dat.label))));

switch MODEFL
    case 'classification'
        OOCVres.predictions = cell(inp1.nclass, MULTILABEL.dim);
        if OOCV.groupmode > 1
            OOCVres.multi_predictions = cell(MULTILABEL.dim,1);
        end
    case 'regression'
        OOCVres.predictions = cell(MULTILABEL.dim,1);
end

if isfield(analysis,'rootdir') && exist(analysis.rootdir,'dir')
    inp1.rootdir = fullfile(analysis.rootdir,analysis.params.TrainParam.SVM.prog,inp1.oocvname);
else
    inp1.rootdir = fullfile(pwd,analysis.params.TrainParam.SVM.prog,inp1.oocvname);
end

if ~exist(inp1.rootdir,'dir'), mkdir(inp1.rootdir); end

% Loop through modalities
for i = 1:inp1.nF
    
    % **************************** ANALYSIS SETUP *****************************
    inp2 = nk_SetFusionMode2(dat, analysis, F, nF, i, inp1.oocvind);
    inp = catstruct(inp1,inp2);
    inp.loadGD = true;
    
    for j = 1:MULTILABEL.dim
	
        strOOCVfile = fullfile(inp.rootdir,[stranalysis inp.varstr '_lb' num2str(j) '_OOCVresults-' num2str(inp.oocvind) '_ID' dat.id '.mat']);
		if exist(strOOCVfile,'file') && inp.ovrwrt==2
			fprintf('\nLoading:\n%s',strOOCVfile);
			load(strOOCVfile)
		else
			if MULTILABEL.flag && MULTILABEL.dim>1
				fprintf('\n\n');cprintf('*black','====== Working on label #%g ====== ',j);
				inp.curlabel = j;
			else
				inp.curlabel = 1;
			end
			if strcmp(MODEFL,'classification')
				switch OOCV.groupmode
					case 1
						inp.multiflag = 0;
						[ijOOCV.BinResults, ijOOCV.FileNames, ijOOCV.RootPath] = nk_OOCV(inp);
					case 2
						inp.multiflag = 1;
						[ijOOCV.MultiResults, ijOOCV.FileNames, ijOOCV.RootPath] = nk_OOCV(inp);
					case 3
						inp.multiflag = 0;
						[ijOOCV.BinResults, ijOOCV.FileNames, ijOOCV.RootPath] = nk_OOCV(inp);
						inp.multiflag = 1;
						[ijOOCV.MultiResults, ijOOCV.FileNames, ijOOCV.RootPath] = nk_OOCV(inp);
				end
            else
                inp.multiflag = 0;
				[ijOOCV.RegrResults, ijOOCV.FileNames, ijOOCV.RootPath] = nk_OOCV(inp);
			end
			fprintf('\nSaving:\n%s',strOOCVfile);
			save(strOOCVfile,'ijOOCV', 'OOCV');	
        end
        
        switch MODEFL
            case 'classification'
                if isfield(ijOOCV,'BinResults')
                    for curclass=1:inp.nclass
                        OOCVres.predictions{curclass, j} = [ OOCVres.predictions{curclass, j} ijOOCV.BinResults.BinCV2Predictions_DecisionValues{curclass} ];
                    end
                    OOCVres.BinResults{j} = ijOOCV.BinResults;
                end
                if isfield(ijOOCV,'MultiResults') 
                    for curclass=1:inp.nclass
                        OOCVres.predictions{curclass, j} = [ OOCVres.predictions{curclass, j} ijOOCV.MultiResults.BinCV2Predictions_DecisionValues{curclass}];
                    end
                    OOCVres.multi_predictions{j} = [ OOCVres.multi_predictions{j} ijOOCV.MultiResults.MultiCV2PredictionsLL ];
                    OOCVres.MultiResults{j} = ijOOCV.MultiResults;
                end
                
            case 'regression'
                OOCVres.predictions{j} = [ OOCVres.predictions{j} ijOOCV.RegrResults.CV2PredictedValues ];
                OOCVres.RegrResults{j} = ijOOCV.RegrResults;
        end
        OOCVres.FileNames{i,j} = ijOOCV.FileNames;
    end
    OOCVres.RootPath{i} = ijOOCV.RootPath;
end
clear inp1 inp2
OOCVres =  nk_OOCVMeta(OOCVres, inp);
hu = findobj('Tag','OOCV');
if ~isempty(hu), delete(hu); end
