function [act, inp] = nk_MLOptimizerPrep(act, inp, parentstr)
% =========================================================================
% [act, inp] = nk_MLOptimizerPrep(act, inp, parentstr)
% =========================================================================
% This function allows the interactive and batch use of the ML training and 
% cross-validation module of NM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 10/2018

global CV NM

na_str = '?';
if ~exist('inp','var') || isempty(inp),
    inp = struct('analind',1, ...
                    'lfl',1, ...
                    'preprocmat',[], ...
                    'gdmat',[], ...
                    'gdanalmat', [], ...
                    'varstr', [], ...
                    'concatfl', [], ...
                    'ovrwrt', 2, ...
                    'update', true, ...
                    'batchflag', false);
     inp = nk_GetAnalModalInfo_config(NM, inp);
     [ ix, jx ] = size(CV(1).TrainInd);
     inp.GridAct = false(ix,jx);
end

%% Configure menu
if isfield(inp,'analind'), 
    if numel(inp.analind)<2,
        AnalSelStr = sprintf('Analysis %g', inp.analind);   
    else
        AnalSelStr = sprintf('%g Analyses: %s',numel(inp.analind), strjoin(cellstr(num2str(inp.analind'))',', '));
    end
else, 
    AnalSelStr = na_str;
end
  
if numel(NM.analysis)>1
    AnalSelectStr = sprintf('Choose analyses to work on [ %s ]|', AnalSelStr);  
    AnalSelectAct = 1;
else
    AnalSelectStr = []; AnalSelectAct = [];
end
% Initialize analysis with current training parameters
analysis = NM.analysis{inp.analind};
if ~exist('inp','var') || isempty(inp)
    if isfield(analysis,'GDdims')
        fprintf('\n')
        cprintf('*red','******************************************************* \n')
        cprintf('*red','ATTENTION: \n')
        cprintf('*red','Previous analysis detected \n')
        cprintf('*red','If you proceed this analysis will be deleted. \n')
        cprintf('*red','******************************************************* \n')
        procflag = nk_input('Proceed anyway',0,'yes|no',[1,0], 2);
        if ~procflag, return; end
        analysis = rmfield(analysis,'GDdims');
        fprintf('\n');cprintf('*red','Training and CV results deleted!'); 
        if isfield(analysis,'visdata'), 
            analysis = rmfield(analysis,'visdata'); 
            fprintf('\n');cprintf('*red','Visualization results deleted!'); 
        end
        if isfield(analysis,'OOCV'),
            analysis = rmfield(analysis,'OOCV'); 
            fprintf('\n');cprintf('*red','Independent test data results deleted!'); 
        end
        analysis.status = 0;
    end
    
end

if ~isempty(analysis)
    
    %% Prepare some useful filename and fieldname definitions
    inp.tdatatypes = {'CVdatamat','CVdatamat','CVresults','CVdimanalysis'};
    inp.sdatatypes = {'','PreprocData','CVdatamat','CVresults'};          
    inp.sfieldnames = {'','preprocmat','gdmat','gdanalmat'};    
    
    if isempty(CV)
         nk_SetupGlobVars2(NM.analysis{inp.analind(1)}.params, 'setup_main', 0); 
    end
    [ix, jx] = size(CV(1).TrainInd); 
    
    
    if isempty(inp.preprocmat),     inp.preprocmat = cell(inp.nF,1); end
    if isempty(inp.gdanalmat),    inp.gdanalmat = cell(inp.nF,1); end
    
    %% Compute from scratch or use pre-computed datamats from different stages of the analysis process?
    
    LFL_opts  = {'Compute from scratch (no CVdatamats or Preprocdata-MATs available)', ...
                    'Compute using existing Preprocdata-MATs (no CVdatamats available)', ...
                    'Assemble analysis using existing CVdatamats', ...
                    'Assemble analysis using existing CVresult-MATs'};
    
    lflcnt = 0;
    if numel(inp.varind)>1,
        LFL_opts(1)=[]; if inp.lfl == 1; inp.lfl = 2; end
        if inp.lfl > 1, lflcnt = 1; end
    end
                
    ModeStr   = sprintf('Operation mode of ML training module [ %s ]|',LFL_opts{inp.lfl-lflcnt});          ModeAct = 2;

    if inp.lfl>1
        % precomputed
        nTargFiles = na_str; mT=1;
        if isfield(inp,inp.sfieldnames{inp.lfl}) && ~isempty(inp.(inp.sfieldnames{inp.lfl})), 
            if iscell(inp.(inp.sfieldnames{inp.lfl}){1})
                nT = sum(sum(~cellfun(@isempty, inp.(inp.sfieldnames{inp.lfl}){1})));
                if numel(inp.(inp.sfieldnames{inp.lfl})) == numel(inp.varind)
                    mT=numel(inp.(inp.sfieldnames{inp.lfl}));
                elseif nT>0 && ~iscell(inp.(inp.sfieldnames{inp.lfl}){1}{1})
                    mT = size(inp.(inp.sfieldnames{inp.lfl}){1}{1},1); 
                else
                    mT = size(inp.(inp.sfieldnames{inp.lfl}),1); 
                end
            else
                nT = sum(sum(~cellfun(@isempty, inp.(inp.sfieldnames{inp.lfl}))));
                if nT>0, 
                    if iscell(inp.(inp.sfieldnames{inp.lfl}))
                        mT = 0;
                    else
                        mT = size(inp.(inp.sfieldnames{inp.lfl}{1}),1); 
                    end
                end
            end
            if mT>1 && nT>0, nvarstr = sprintf('%gx',mT); else, nvarstr = ''; end
            nTargFiles = sprintf('%s%g selected', nvarstr, nT); 
        end     
        FilesStr = sprintf('Specify %s files [ %s ]|', inp.sdatatypes{inp.lfl}, nTargFiles);        FilesAct = 3;
    else
        FilesStr = ''; FilesAct = [];
    end
    % from scratch
    OVRWRT_opts = { sprintf('Overwrite %s files',inp.tdatatypes{inp.lfl}), ...
                sprintf('Do not overwrite %s files',inp.tdatatypes{inp.lfl})};      

    OverWriteStr = sprintf('Overwrite existing %s files [ %s ]|', ...
        inp.tdatatypes{inp.lfl}, OVRWRT_opts{inp.ovrwrt}) ;                                         OverWriteAct = 4; 
    
    %% Retrieve CV2 partitions to operate on
    if ~isfield(inp,'GridAct'), inp.GridAct = false(ix,jx); end;                                              
    GridSelectStr = sprintf('Select CV2 partitions to operate on [ %g selected ]|', ...
                                                                            sum(inp.GridAct(:)));  GridSelectAct = 5;
                                                                        
    %% Build interactive menu
    menustr = [ AnalSelectStr ...
                ModeStr ...
                FilesStr ...
                OverWriteStr ...
                GridSelectStr ];

    menuact = [ AnalSelectAct ...
                ModeAct ...
                FilesAct ...
                OverWriteAct ...
                GridSelectAct ];       

    disallow = false;

    %% Check whether all parameters are available
    if ~sum(inp.GridAct(:)) || isempty(inp.analind), disallow = true; end
    if inp.lfl>1, if ~isfield(inp, inp.sfieldnames{inp.lfl}) || isempty(inp.(inp.sfieldnames{inp.lfl})), disallow = true; end; end

    if ~disallow, menustr = [menustr 'PROCEED >>>']; menuact = [menuact 6]; end

	% Set run-time flags
    if inp.ovrwrt == 2, ovrwrt = false; else, ovrwrt = true; end
    switch inp.lfl
        case 1
            inp.ovrwrtGD        = ovrwrt; 
            inp.updGD           = true; 
            inp.ovrwrtRes       = true; 
            inp.updRes          = true;
        case 2
            inp.ovrwrtGD        = ovrwrt;
            inp.updGD           = inp.update; 
            inp.ovrwrtRes       = true; 
            inp.updRes          = true;
        case 3
            inp.ovrwrtGD        = false;
            inp.updGD           = inp.update;
            inp.ovrwrtRes       = ovrwrt; 
            inp.updRes          = true; 
        case 4
            inp.ovrwrtGD        = false;
            inp.updGD           = false;
            inp.updRes          = inp.update;
            inp.ovrwrtRes       = false;
    end
	
    %% Display menu and act on user selections
    nk_PrintLogo
    mestr = 'ML Training module run-time configuration'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr);
	if ~inp.batchflag, act = nk_input(mestr, 0, 'mq', menustr, menuact); end
    algostr = getAlgoStr(analysis);
    switch act 
        case 0
            return
        case 1
            if isfield(inp,'analind'), analind = inp.analind; end
            brief = 1;t_act = 1; while t_act>0, [t_act, analind, ~, ~, brief] = nk_SelectAnalysis(NM, 0, 'MAIN INTERFACE >> TRAIN ML MODELS', analind, [], 0, 0, brief); end;
            if ~isempty(analind)
                inp.analind             = analind; 
                inp.preprocmat          = cell(inp.nF,1); 
                inp.gdmat               = [];
                inp.gdanalmat           = []; 
                inp = nk_GetAnalModalInfo_config(NM, inp);   
                [ ix, jx ] = size(CV(1).TrainInd);
                inp.GridAct = false(ix,jx);
            end 
        case 2
            if numel(inp.varind)>1, selLFLopts = 2:4; lfldef = inp.lfl - 1; else, selLFLopts = 1:4; lfldef = inp.lfl; end
            lfl = nk_input('Define run-time mode of ML training module',0,'mq',strjoin(LFL_opts, '|'), selLFLopts, lfldef);
            if numel(inp.varind)>1 && lfl == lfl+1; end
            if lfl, inp.lfl = lfl; end        
        case 3
			%% Prepare for runtime:
			% This section of the code sets the options and file masters for the ML optimization procedure
			% Basically we distinguish between:
            switch inp.lfl 
                case 1 % Running from scratch
                    inp.preprocmat          = cell(inp.nF,1); 
                    inp.gdmat               = [];
                    inp.gdanalmat           = []; 
                    inp.GridAct             = true(ix,jx); 
                case 2 % Running with precomputed Preprocdatamats
                    [preprocmat, emptfl]    = nk_GenPreprocMaster2(NM.id, CV(1), [],  NM.analysis{inp.analind}.rootdir, [], [], inp.varind, inp.varstr, inp.concatfl);
                    if ~emptfl && ~isempty(preprocmat), 
                        inp.preprocmat      = preprocmat;
                        inp.gdmat           = [];
                        inp.gdanalmat       = []; 
                        inp.GridAct         = ~cellfun(@isempty,preprocmat{1}); 
                    end
                case 3 % Running with precomputed CVdatamats (aggregation run --- allowing to tweak some post-training model selection options)
                   
                    [gdmat, emptfl]         = nk_GenCVdataMaster2(NM.id, CV(1), [], fullfile(NM.analysis{inp.analind}.rootdir, algostr), [], [], inp.varind, inp.varstr, inp.concatfl);
                    if ~emptfl && ~isempty(gdmat), 
                        inp.preprocmat      = cell(inp.nF,1);
                        inp.gdmat           = gdmat; 
                        inp.GridAct         = ~cellfun(@isempty,gdmat{1}); 
                    end 
                case 4 % Running with existing CVresults (simple aggregation run --- basically for multiple modalities)
                    
                    [gdanalmat, emptfl]     = nk_GenCVresultsMaster(NM.id,[], fullfile(NM.analysis{inp.analind}.rootdir, algostr)); 
                    if ~emptfl && ~isempty(gdanalmat), 
                        inp.preprocmat      = cell(inp.nF,1);
                        inp.gdmat           = [];
                        inp.gdanalmat       = gdanalmat;
                        load(gdanalmat{1}); 
                        inp.GridAct         = GDanalysis.GridAct;
                        clear GDanalysis
                    end
            end
        case 4
            if inp.ovrwrt ==1, inp.ovrwrt=2; elseif inp.ovrwrt == 2, inp.ovrwrt = 1; end
        case 5
            t_act = 1; while t_act > 0 && t_act < 11, [ t_act, inp.GridAct ] = nk_CVGridSelector(ix,jx, inp.GridAct, 0); end
        case {6,7}
            if inp.lfl ~=2, inp.preprocmat = []; end
			% act==7 is the automation option. Make sure that inp is properly defined 
			% To this end create an inp structure based on 
            nA = 1; if numel(inp.analind)>1, nA = numel(inp.analind); end
            for i=1:nA
                nk_SetupGlobVars2(NM.analysis{inp.analind(i)}.params, 'setup_main', 0); 
                NM.runtime.curanal = inp.analind(i);
                inp.analysis_id = NM.analysis{inp.analind(i)}.id;
                NM.analysis{inp.analind(i)} = MLOptimizerPrep(NM, NM.analysis{inp.analind(i)}, inp);
                nk_SetupGlobVars2(NM.analysis{inp.analind(i)}.params, 'clear', 0); 
            end
            h = findobj('Tag','PrintCVBarsBin'); if ~isempty(h), delete(h); end
    end
end

function analysis = MLOptimizerPrep(dat, analysis, inp1)

global PARMODE MULTI SAV MODEFL CV PREPROC RAND FUSION META SVM

nk_SetupGlobVars2(analysis.params, 'setup_strat', 0, inp1.varind(1)); 
strout = nk_Preprocess_StrCfg([], []);

% Define # of classifiers to train (1 for multi-group classification & regression)
inp1.nclass = 1; if strcmp(MODEFL,'classification') && RAND.Decompose ~= 9, inp1.nclass = length(CV(1).class{1,1}); end
inp1.probflag = false; 

% **************************** ANALYSIS SETUP *****************************
ld = 1; if FUSION.flag == 3, ld = numel(inp1.F); end
inp1.ngroups = max(dat.label);
hx = size(dat.label,2);
analysis.Time                               = zeros(ld,1);
analysis.TrainPerformanceBin                = zeros(ld,inp1.nclass,hx);
analysis.TestPerformanceBin                 = zeros(ld,inp1.nclass,hx);
analysis.TestPerformanceBinPermAggr         = zeros(ld,inp1.nclass,hx);
analysis.TestPerformanceBinPermAggrMajVote  = zeros(ld,inp1.nclass,hx);
analysis.TrainPerformanceMean               = zeros(ld,1);
analysis.TestPerformanceMean                = zeros(ld,1);
analysis.TestPerformanceMeanPermAggr        = zeros(ld,1);

if ~isempty(MULTI) && MULTI.flag && inp1.nclass > 1
    analysis.TrainPerformanceMulti          = zeros(ld,1);
    analysis.TestPerformanceMulti           = zeros(ld,1);
    analysis.TestPerformanceMultiPermAggr   = zeros(ld,1);
    analysis.TestPerformanceMultiPermAggrGroup = zeros(ld,inp1.ngroups,1);
    inp1.multiflag = true;
else
    inp1.multiflag = false; 
end

PARMODE = 0; 
analysis.status = 0;
GDdims = cell(inp1.nF,1);
preprocmat = inp1.preprocmat;

switch FUSION.flag
    case 3
        mthstr = 'DECFUSE';
    otherwise
        mthstr = SVM.prog;
end

% Create analysis directory if needed
if ~isfield(inp1,'rootdir') || isempty(inp1.rootdir) || ~exist(inp1.rootdir,'dir')==7
    inp1.rootdir = fullfile(pwd,mthstr);
else
    inp1.rootdir = fullfile(analysis.rootdir,mthstr); 
end

if ~exist(inp1.rootdir,'dir'), mkdir(inp1.rootdir);end
 
for i = 1:inp1.nF
    
    %% Get Training / CV data (Y) & Build modality suffix
    inp2  = nk_SetFusionMode2(dat, analysis, inp1.F, inp1.nF, i);
    inp   = catstruct(inp1,inp2); clear inp2
    if isfield(dat,'time'), inp.time2event = dat.time; end
    if inp.lfl == 2, inp.preprocmat = preprocmat{i,:,:}; end
    
    inp.stranalysis = SAV.matname; strGDdimsfile = fullfile(inp.rootdir, ...
                        [inp.stranalysis '_CVdimanalysis' strout '_ID' dat.id '.mat']);
    
    PX = nk_ReturnParamChain(PREPROC, 1); PreML = [];
    
    if ~isempty(PX)
        nP = numel(PX);
        if nP>1,
            Params = []; Params_desc=[]; ModalityVec = [];
            for n=1:nP
                Params = [Params PX(n).Params];
                for m=1:numel(PX(n).Params_desc)
                    PX(n).Params_desc{m} = sprintf('%s_{M%g}',PX(n).Params_desc{m},n);
                end
                Params_desc = [Params_desc PX(n).Params_desc] ;
                ModalityVec = [ModalityVec n*ones(1,numel(PX(n).Params))];
            end
            PreML.Params_desc   = Params_desc;
            PreML.Params        = Params;
            PreML.ModalityVec   = ModalityVec;
        else
            PreML.Params_desc   = PX.Params_desc;
            PreML.Params        = PX.Params;
            PreML.ModalityVec   = ones(1,numel(PX.Params));
        end
    end

    maxacc = zeros(inp.nclass, 1); maxtestacc = zeros(inp.nclass,1 );
    
    % *********************************** GO **********************************
    if ~isempty(inp.gdanalmat{i})
        % Load precomputed CVresults file
        load(inp.gdanalmat{i}); GDdims{i} = GDanalysis;
    else
        % Train classifier / predictor at current dimensionality
        GDdims{i} = nk_MLOptimizer_main(inp, dat.id, PreML);
    end
    if ~inp.batchflag
        GDdims{i}.datadescriptor = dat.datadescriptor{inp.tF};
        % Retrieve best accuracies of binary classifiers from grid structure
        for j = 1:inp.nclass
            maxacc(j) = GDdims{i}.best_CVperf{j};
            analysis.TrainPerformanceBin(i,j) = maxacc(j);
            maxtestacc(j) = GDdims{i}.best_TSperf{j};
            analysis.TestPerformanceBin(i,j) = maxtestacc(j);
            switch MODEFL
                case 'regression'
                    analysis.TestPerformanceBinPermAggr(i,j,:)        = GDdims{i}.Regr.costfun_crit;
                case 'classification'
                    analysis.TestPerformanceBinPermAggr(i,j,:)        = GDdims{i}.BinClass{j}.costfun_crit;
                    if RAND.Decompose ~= 9
                        analysis.TestPerformanceBinPermAggrMajVote(i,j,:) = GDdims{i}.BinClass{j}.prob_contigency.BAC;             
                    end
            end
        end

        % 1) Calculate maximum performance for current dimensionality by means of 
        % averaging the binary classifiers 
        meanmaxacc                              = mean(maxacc);
        analysis.TrainPerformanceMean(i)        = meanmaxacc;
        analysis.TestPerformanceMean(i)         = mean(maxtestacc(maxtestacc~=0));
        analysis.TestPerformanceMeanPermAggr(i) = mean(analysis.TestPerformanceBinPermAggr(i,:));
        if ~isempty(MULTI) && MULTI.flag && inp.nclass > 1
            analysis.TrainPerformanceMulti(i) = GDdims{i}.best_MultiCVperf;
            analysis.TestPerformanceMulti(i) = GDdims{i}.best_MultiTSperf;
            analysis.TestPerformanceMultiPermAggr(i) = GDdims{i}.MultiClass.accuracy;
            for z=1:inp.ngroups
                analysis.TestPerformanceMultiPermAggrGroup(i,z) = GDdims{i}.MultiClass.class{z}.BAC;
            end
        end
    end
end

% Store analysis results in analysis structure and set status to 1
analysis.GDdims = GDdims;
if ~isempty(META) && META.flag, 
    IN.labels = inp.labels;
    IN.nclass = inp.nclass;
    IN.ngroups = unique(inp.labels); IN.ngroups(isnan(IN.ngroups))=[]; IN.ngroups = numel(IN.ngroups);
    METAres = nk_MLOptimizerMeta( analysis, IN ); 
    if ~isempty(METAres), analysis.META = METAres; end
end
    
analysis.status = 1;
% Save analysis to disk
if ~inp.batchflag
    fprintf('\nSaving %s', strGDdimsfile);
    save(strGDdimsfile,'analysis');
end

function [DR, BINMOD, FEATSEL, CLUST, COVAR] = get_preproc_params(PREPROC)

DR=[]; FEATSEL=[]; CLUST=[]; COVAR=[]; BINMOD=[]; 
if isfield(PREPROC,'FEATSEL'), FEATSEL = PREPROC.FEATSEL; end
if isfield(PREPROC,'BINMOD'), BINMOD = PREPROC.BINMOD; end
if isfield(PREPROC,'ACTPARAM')
    for i=1:numel(PREPROC.ACTPARAM) 
       if isfield(PREPROC.ACTPARAM{i},'DR')
           DR = PREPROC.ACTPARAM{i}.DR;
       elseif isfield(PREPROC.ACTPARAM{i},'CLUST')
           CLUST = PREPROC.ACTPARAM{i}.CLUST;
       elseif isfield(PREPROC.ACTPARAM{i},'COVAR')
           COVAR = PREPROC.ACTPARAM{i}.COVAR;  
       end
    end
end