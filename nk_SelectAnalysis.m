function [act, analdim, dat, showmodalvec, brief] = nk_SelectAnalysis(dat, newflag, ...
    titlestr, analind, modflag, complflag, showmodalvec, brief)
% =========================================================================
% [act, analdim, dim, showmodalvec] = nk_SelectAnalysis(dat, newflag, ...
%                     titlestr, analind, modflag, complflag, showmodalflag)
% =========================================================================
% Helper function to select an analysis for classification, permutation,
% visualization ...
%
% Inputs:
% -------
% dat               : NeuroMiner structure
% newflag           : Generate new analysis
% titlestr          : Action string to show in the multi selection menu
% analind           : Start with a given analysis (default = 1)
% modflag           : Overwrite Training parameters of dat with current
%                     analysis params
% complflag         : Show only analyses that have been computed
% showmodalvec      : Specific list of modalities to show to user  
%
% Outputs:
% --------
% act               : Selected menu item
% analind           : Selected analysis index
% dat               : (modified) NM structure
% showmodalvec      : (modified) list of modality to show to user  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2017

if ~exist('newflag','var')  || isempty(newflag),     newflag = false; end
if ~exist('modflag','var')  || isempty(modflag),     modflag = false; end
if ~exist('complflag','var')|| isempty(complflag),   complflag = false; end
if ~exist('brief','var')    || isempty(brief),       brief = true; end
if ~exist('titlestr','var') || isempty(titlestr),    titlestr='Select'; end

if iscell(dat.Y), nvar = length(dat.Y); else nvar = 1; end

analdim = []; act = 0;

% Select analysis 
if isfield(dat,'analysis')
    
    if complflag 
        complstr = 'COMPLETED ';
        analstatus = nk_GetAnalysisStatus(dat);
        analyses = dat.analysis(analstatus.completed_analyses);
    else
        complstr = '';
        analyses = dat.analysis;
    end
    
    n = length(analyses);
    if brief
        analsel = print_analysis_quickselector(analyses);
        if ~analsel 
            return; 
        elseif analsel > 0
            if complflag
                analselcompl = find(analstatus.completed_analyses);
                analdim = analselcompl(analsel);
                return
            else
                analdim = analsel;
                return
            end
        else
            brief = 0;
        end
    end
    if ~exist('analind','var') || isempty(analind)
        analind=1; 
    elseif analind > n
        analind = n;
    elseif analind < 1
        analind = 1;
    end
    
    clc
    cprintf('*red','%s \n',repmat('*',1,size(titlestr,2)));
    cprintf('*red','%s \n', titlestr)
    cprintf('*red','%s \n',repmat('*',1,size(titlestr,2)));
    fprintf('\n')
    if n == 1
        analstr = sprintf('1 %sANALYSIS',complstr);
    else
        analstr = sprintf('CURRENT ANALYSIS: %g (of %g %sANALYSES)', analind, n, complstr);
    end
    actstr = ['ANALYSIS MANAGER: ' analstr];
    nF=1;
    % Retrieve some information
    if isfield(analyses{analind},'params') && ~isempty(analyses{analind}.params)
        
        skipfl = false; showmodalmax = 3;
        showmodalvec = 1:nvar;
        params = analyses{analind}.params;
        FUSION = params.TrainParam.FUSION;
        if isempty(FUSION), FUSION.flag = 0; end
        if nvar > 1
            if params.TrainParam.STACKING.flag == 1
                cprintf('*black', 'STACKING ANALYSIS OPERATES ON ANALYSES %s ', strjoin(string(params.TrainParam.STACKING.sel_anal),', ')); 
                fprintf('\n\tStacking descriptors: %s', strjoin(params.TrainParam.STACKING.featnames,','))
                predstr = {'CV1 training data','CV1 test data'};
                fprintf('\n\tPrediction extraction for stacking from %s ', predstr{params.TrainParam.STACKING.mode});
            else
                nF = numel(FUSION.M);
                if nF>1,
                    strmod = sprintf('%g MODALITIES',nF);
                else
                    strmod = 'MODALITY';
                end
                cprintf('*black', 'ANALYSIS OPERATES ON %s ', strmod); 
                for j=1:nF, 
                    if j > showmodalmax, fprintf(' ...'); break; end;
                    cprintf('*black','#%g ', FUSION.M(j)); 
                end
            end
            fprintf('\n') 
        end
        e = nk_GetParamDescription2(dat, params,'cv');
        %nV = numel(showmodalvec);
        
        % Loop through variates
        for jj=1:nF
            
            j = showmodalvec(jj);

            if FUSION.flag == 1 && j>1
                continue; 
            elseif FUSION.flag == 3
                params.TrainParam = analyses{analind}.params.TrainParam.STRAT{j};
            else
                params.TrainParam = analyses{analind}.params.TrainParam;
            end
            
            switch FUSION.flag
                case {0,2}
                    % Get info about modality j
                    d = nk_GetParamDescription2(dat, params, 'VarDesc', [], FUSION.M(j));
                    d = nk_GetParamDescription2(dat, params.TrainParam.PREPROC{FUSION.M(j)},'PreProc', d, j);
                case 1
                    % Get info about all modalities
                    d = nk_GetParamDescription2(dat, params,'VarDesc', [], FUSION.M);
                    d = nk_GetParamDescription2(dat, params.TrainParam.PREPROC{j},'PreProc', d, 1);
                case 3
                    % Get info about modality j
                    d = nk_GetParamDescription2(dat, params,'VarDesc', [], FUSION.M);
                    d = nk_GetParamDescription2(dat, params.TrainParam.PREPROC,'PreProc', d);                    
            end
            
            if FUSION.flag == 1
                mxl = 0;
                for i=1:numel(d.datadescriptor)
                    mdlstr = sprintf('MODALITY %g : %s \n', FUSION.M(i), d.datadescriptor{i});
                    cprintf('*blue',mdlstr);
                    mxli = size(mdlstr,2);
                    if mxli > mxl, mxl = mxli; end
                end
            else
                if numel(d.datadescriptor) > 1,
                    datdesc = d.datadescriptor{j};
                else
                    datdesc = d.datadescriptor{1};
                end
                mdlstr = sprintf('MODALITY %g : %s \n', FUSION.M(j), datdesc);
                cprintf('*blue',mdlstr);
                mxl = size(mdlstr,2);
            end
            
            cprintf('*blue','%s \n',repmat('*',1,mxl));
            cprintf('*black','Preprocessing: \n'); 
            
            if strcmp(dat.modeflag,'classification')
                fprintf('\t* %s \n', d.PREPROC.groupmode);
            else
                fprintf('\t* %s \n', d.PREPROC.targetscaling);
            end
            
            for k=1:numel(d.PREPROC.preprocact)
                fprintf('\t* Step %g: %s \n', k, d.PREPROC.preprocact{k}); 
            end 
            
            if FUSION.flag == 3, print_modalitydata(dat, params, 1); end
            
            if jj > showmodalmax, 
                cprintf('*blue','%s \n',repmat('-',1,mxl)); 
                cprintf('*red','>>> %g further modalities included in this analysis... \n', nF-j+1);
                skipfl = true;
                break; 
            else
                if j<nvar, cprintf('*blue','%s \n',repmat('-',1,mxl)); end
            end
        end
    end
    if ~skipfl && FUSION.flag ~= 3, print_modalitydata(dat, params, FUSION.M); end
    if nF > 1, cprintf('*blue','%s \n',repmat('=',1,100)); end
    cprintf('*black','Cross-Validation: '); fprintf('\n\t* %s\n\n', e.cv);
    cprintf('*black','Results: ');
    istr = [];
    
    if analyses{analind}.status
        ClassHd         = 'Predictor/Classifier        ';
        ClassBl         = blanks(length(ClassHd));
        if size(analyses{analind}.TrainPerformanceBin,1)>1
            latfus = [' (+' num2str(nvar-1) ' M.)']; 
        else
            latfus = '';
        end
        TrBinStrHd      = 'Mean CV1';
        TsBinStrHd      = 'Mean CV2';
        TsBinAggrStrHd  = 'OOT-CV2';
        hdrstr = sprintf('\n\n\t%s\t%s%s\t%s%s\t%s%s\n', ClassHd, TrBinStrHd, latfus, TsBinStrHd, latfus, TsBinAggrStrHd, latfus);
        fprintf('%s',hdrstr);
        fprintf('\t%s\n',repmat('=',1,length(hdrstr)+4));
        
        for q=1:size(analyses{analind}.TrainPerformanceBin,2)
            TrBinStr     = analyses{analind}.TrainPerformanceBin(1,q);
            TsBinStr     = analyses{analind}.TestPerformanceBin(1,q);
            TsBinAggrStr = analyses{analind}.TestPerformanceBinPermAggr(1,q);
            switch dat.modeflag
                case 'classification'
                    ClassStr    = analyses{analind}.params.cv.class{1,1}{q}.groupdesc;
                case 'regression'
                    ClassStr    = 'Regression';
            end
            ClassStrL    = length(ClassStr);
            ClassHd      = ClassBl; ClassHd(1:ClassStrL) = ClassStr;
            fprintf('\t%s\t%7.1f%s\t%10.1f%s\t%10.1f%s\n', ClassHd, TrBinStr, latfus, TsBinStr, latfus, TsBinAggrStr, latfus);
        end
        fprintf('\t%s\n',repmat('-',1,length(hdrstr)+4));
        istr=[];
        if params.TrainParam.MULTI.flag

            MultiCV1ParamStr = num2str(analyses{analind}.TrainPerformanceMulti,'%1.1f');
            MultiCV2ParamStr = num2str(analyses{analind}.TestPerformanceMulti,'%1.1f');
            MultiCV2ParamOOTStr = num2str(analyses{analind}.TestPerformanceMultiPermAggr,'%1.1f');
            ClassStr     = 'Multi-Group: Overall';
            ClassStrL    = length(ClassStr);
            ClassHd      = ClassBl; ClassHd(1:ClassStrL) = ClassStr;
            if size(MultiCV1ParamStr,1)>1
                latfus = [' (+' num2str(nvar-1) ' M.)']; 
                istr = [istr ...
                    '\t\t\t\t' ClassHd '  ' ...
                       MultiCV1ParamStr(1,:) latfus '     ' ...
                       MultiCV2ParamStr(1,:) latfus '     ' ...
                       MultiCV2ParamOOTStr(1,:) latfus '\n'];
            else
                latfus = '';
                istr = [istr ...
                    '\t\t\t\t' ClassHd '  ' ...
                       MultiCV1ParamStr '      ' ...
                       MultiCV2ParamStr '      ' ...
                       MultiCV2ParamOOTStr '\n'];
            end
            for q=1:size(analyses{analind}.TestPerformanceMultiPermAggrGroup,2)
                MTsMultiAggrGrpStr = [num2str(analyses{analind}.TestPerformanceMultiPermAggrGroup(1,q),'%1.1f') latfus];
                MClassStr     = ['Multi-Group: ' dat.groupnames{q}];
                MClassStrL    = length(MClassStr);
                MClassHd      = ClassBl; MClassHd(1:MClassStrL) = MClassStr;
                istr = [istr ...
                    '\t\t\t\t' MClassHd '  '...
                    '          ' ...
                    '          ' ... 
                    MTsMultiAggrGrpStr '\n'];
            end
        end

    else
        istr = [istr ... 
            'empty\n'];
    end
    
    if skipfl
       menustr = '|Show me other modalities ...';
       menudat = 5;
    else
       menustr = [];
       menudat = [];
    end
    
    menustr = [menustr '|Select current analysis'];
    menudat = [menudat 4];
    
    if newflag
        menustr = [menustr '|Generate new analysis'];
        menudat = [menudat 3];
    end
    
    if n > 1
        if modflag 
            menustr = [menustr '|Overwrite current parameter template with parameters of this analysis'];
            menudat = [menudat 6];
        end
        if analind == n && n ~= 1
            menustr = [menustr '|<= Go to previous analysis'];
            menudat = [menudat 1];
            if n > 2 && analind > 2,
                menustr = [menustr '|<<= Go to first analysis'];
                menudat = [menudat 7];
            end
        elseif analind == 1 && n ~= 1
            menustr = [menustr '|Go to next analysis =>'];
            menudat = [menudat 2];
            if n > 2 && analind < (n - 1),
                menustr = [menustr '|Go to last analysis =>>'];
                menudat = [menudat 8];
            end
        else
            menustr = [menustr '|Go to next analysis =>|<=Go to previous analysis'];
            menudat = [menudat 2,1];
            if n > 2 && ( analind > 2 && analind < (n - 1) )
                if analind > 2,
                    menustr = [menustr '|<<= Go to first analysis'];
                    menudat = [menudat 7];
                end
                if analind < (n - 1),
                    menustr = [menustr '|Go to last analysis =>>'];
                    menudat = [menudat 8];
                end
            end
        end
        if n>3
             menustr = [menustr '|Jump to analysis ...'];
             menudat = [menudat 9];
        end
    end
 
    act = nk_input(actstr,0,'mq',menustr,menudat);
    
    switch act
        case 1
            analdim = analind - 1;
        case 2
            analdim = analind + 1;
        case 3
            analdim = n + 1; act = 0;
        case 4
            analdim = analind; act = 0;
        case 5
            showmodalvec = nk_SelectVariateIndex(dat,0);
            analdim = analind;
        case 6
           surefl = nk_input('This will the overwrite the current parameter template! Are your sure?',0,'yes|no',[1,0],2);
           if surefl
               dat.cv = analyses{analind}.params.cv;
               dat.TrainParam = analyses{analind}.params.TrainParam;
               dat.SVM = analyses{analind}.params.SVM;
               if isfield(analyses{analind}.params,'RVM')
                   dat.RVM = analyses{analind}.params.RVM;
               end
               if isfield(analyses{analind}.params,'MKLRVM')
                   dat.MKLRVM = analyses{analind}.params.MKLRVM;
               end
           end
           act = 0;
        case 7 
            analdim = 1;
        case 8
            analdim = n;
        case 9
            t_analdim = nk_input('Jump to analysis',0,'w1');
            if t_analdim > 0 && t_analdim <= numel(analyses), analdim = t_analdim; end
    end
else
    analdim = 1;
end

function print_modalitydata(dat, params, varind)

    e = nk_GetParamDescription2(dat, params.TrainParam.RFE,'FeatFlt');
    e = nk_GetParamDescription2(dat, params.TrainParam.RFE,'FeatWrap',e);
    e = nk_GetParamDescription2(dat, params.TrainParam,'multiclass',e);
    e = nk_GetParamDescription2(dat, params.TrainParam,'GridParam',e);
    e = nk_GetParamDescription2(dat, params.TrainParam,'ParamComb',e, varind);
    e = nk_GetParamDescription2(dat, params.TrainParam,'SVMprog',e);
    e = nk_GetParamDescription2(dat, params.TrainParam,'classifier',e);
    e = nk_GetParamDescription2(dat, params.TrainParam,'kernel',e);
    
    % Generate analysis description
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    if params.TrainParam.RFE.Filter.flag, 
        cprintf('*black','Feature selection (Filter): '); fprintf('\n\t* %s\n',e.FilterMode);
        fprintf('\t* %s\n', e.FilterMethod); 
    end
    if params.TrainParam.RFE.Wrapper.flag,
        cprintf('*black','Feature selection (Wrapper): '); fprintf('\n\t* %s\n',e.WrapperStr);    
        fprintf('\t* %s\n', e.WrapperMethod); 
    end
    if params.TrainParam.MULTI.flag, 
        cprintf('*black','Multi-group classification: '); fprintf('\n\t* %s\n', e.multiclass);
    end
    cprintf('*black','Machine Learning Method: '); fprintf('\n\t* %s, %s, %s\n', e.prog, e.classifier, e.kernel);
    if e.preML_nCombs>0
        cprintf('*black','Preprocessing Optimization: %g parameter combinations', e.preML_nCombs); 
        for i=1:numel(e.preML), fprintf('\n\t(%g) %s', i, e.preML{i}); end; fprintf('\n');
    else
        cprintf('*black','No optimization of preprocessing parameters \n')
    end
    if e.ML_nCombs>0
        cprintf('*black','ML Optimization: %g parameter combinations',e.ML_nCombs); 
        for i=1:numel(e.ML), fprintf('\n\t(%g) %s',i, e.ML{i}); end; fprintf('\n');
    else
        cprintf('*black','No optimization of ML parameters \n')
    end
    
function analind = print_analysis_quickselector(analyses)
    
nk_PrintLogo
fprintf('\n\t'); cprintf('blue*', '============================================= '); 
fprintf('\n\t'); cprintf('black*','***        Quick Analysis Selector        *** ');
fprintf('\n\t'); cprintf('blue*', '============================================= ');
for i=1:numel(analyses)
    fprintf('\n\t** Analysis [ %3g ]: ID : %s', i, analyses{i}.id);
end
fprintf('\n');
analsel = nk_input(sprintf('Type the number of analysis you want to choose,\n\t''S'' for analysis sequence\n\t''M'' for detailed analysis information,\n\tor some other string to return'),0,'sq');
if any(strcmp({'S','s'},analsel))
    analind = nk_input('Type sequence of analyses you want to work on',0,'e');
elseif ~any(strcmp({'M','m'}, analsel)) && ~isnan(str2double(analsel)) 
    analind = str2double(analsel);
    if analind < 1 || analind > numel(analyses), analind = 0; end
elseif any(strcmp({'M','m'}, analsel))
    analind = -1;
else
    analind = 0;
end
        
    
    
