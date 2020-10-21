function [GD, MD, DISP] = nk_GridSearchHelper(GD, MD, DISP, i, nclass, ngroups, CV1PerfData, CV2PerfData, models)

global VERBOSE MULTI CV BATCH SVM MODEFL W2AVAIL RAND MULTILABEL RFE SAV

%%%%%% TRANSFER RESULTS TO GD :

if SAV.savemodel && (exist('models','var') && ~isempty(models))
    MD{i} = models;
end

%% CV1 test data performance
% binary

CV2Param = RFE.CV2Class;
if RFE.Wrapper.flag 
    Param = RFE.Wrapper;
else
    Param = RFE.Filter;
end

if isfield(Param,'EnsembleStrategy')
    Metric = Param.EnsembleStrategy.Metric;
else
    Metric = CV2Param.EnsembleStrategy.Metric;
end

DISP.nclass = nclass; LCOstr = '';

curlabel = MULTILABEL.curdim;

if MULTILABEL.dim > 1
    if RFE.Wrapper.flag
        Perf = CV1PerfData{curlabel}.Wrapper;
    else
        Perf = CV1PerfData{curlabel}.Filter;
    end
    CV2Perf = CV2PerfData{curlabel};
else
    if RFE.Wrapper.flag
        Perf = CV1PerfData.Wrapper;
    else
        Perf = CV1PerfData.Filter;
    end
    CV2Perf = CV2PerfData;
end

switch Metric

    case 1 % Mean of ensembles (Targets)

        GD.TR(i,:,curlabel) = Perf.MeanCVHTperf;

    case 2 % Mean of ensembles (Decision Values)

        GD.TR(i,:,curlabel) = Perf.MeanCVHDperf;

end

if (Param.flag && Param.SubSpaceFlag ) && (isfield(Param,'EnsembleStrategy') && Param.EnsembleStrategy.type ~= 9 )

    GD.M_DivT(i,:,curlabel)    = Perf.Ens_MeanTrDiv;
    GD.SD_DivT(i,:,curlabel)   = Perf.Ens_SDTrDiv;
    GD.M_DivV(i,:,curlabel)    = Perf.Ens_MeanCVDiv;
    GD.SD_DivV(i,:,curlabel)   = Perf.Ens_SDCVDiv;
    flg = true;     
else
    flg = false;
end

% multi-group:
if MULTI.flag

   GD.MultiCV1TrPred{i,curlabel}       = Perf.MultiTrPredictions;
   GD.MultiCV1CVPred{i,curlabel}       = Perf.MultiCVPredictions;
   GD.MultiTR(i,curlabel)              = Perf.MeanMultiCVPerf; 
    switch CV2Param.type 
        case 1 % mean of CV1 ensemble
            GD.MultiTS(i,curlabel) = mean(CV2Perf.MultiCV1Performance(:));
        case 2 % Ensemble of ensembles
            GD.MultiTS(i,curlabel) = CV2Perf.MultiCV2Performance;
    end
    GD.MultiCV2Div(i,curlabel) = CV2Perf.MultiCV2Diversity_Targets;
    GD.MultiCV2DivDec(i,:,curlabel) = CV2Perf.MultiCV2Diversity_DecValues;
    GD.MultiCV2Pred{i,curlabel}   = CV2Perf.MultiCV2Predictions;
    GD.MultiCV2Prob{i,curlabel}   = CV2Perf.MultiCV2Probabilities;
    GD.MultiCV1Pred{i,curlabel}   = nk_cellcat(CV2Perf.MultiCV1Predictions,[],2);
    for g=1:ngroups
        F = repmat({g}, size(CV2Perf.MultiCV1Probabilities));
        GD.MultiCV1Prob{i,g,curlabel} = nk_cellcat(CV2Perf.MultiCV1Probabilities, [], 2, F);
    end
    GD.MultiERR(i,curlabel)    = GD.MultiTR(i) - GD.MultiTS(i);

   if flg
        GD.MultiM_DivT(i,curlabel)     = Perf.Ens_MeanMultiTrDiv;
        GD.MultiSD_DivT(i,curlabel)    = Perf.Ens_SDMultiTrDiv;
        GD.MultiM_DivV(i,curlabel)     = Perf.Ens_MeanMultiCVDiv;
        GD.MultiSD_DivV(i,curlabel)    = Perf.Ens_SDMultiCVDiv;
   end
end

%% CV2 test data performance
% binary:
switch CV2Param.type
    case 1
        % Mean of CV1 ensemble
        GD.TS(i,:,curlabel) = CV2Perf.BinCV1Performance_Mean;
    case 2
        % Ensemble of ensemble decision
        switch Metric                                     
            case 1
                GD.TS(i,:,curlabel) = CV2Perf.binCV2Performance_Targets;
            case 2
                GD.TS(i,:,curlabel) = CV2Perf.binCV2Performance_DecValues;
        end
end

% Decision values / Probabilities on CV2 test data
GD.DS{i,curlabel}          = CV2Perf.binCV1Predictions;
GD.BinPred{i,curlabel}     = CV2Perf.binCV2Predictions_DecValues;
GD.CV2Div(i,:,curlabel)    = CV2Perf.binCV2Diversity_Targets;
GD.CV2DivDec(i,:,curlabel) = CV2Perf.binCV2Diversity_DecValues;

% Model params
GD.FEAT{i,curlabel}        = Perf.SubSpaces;
if flg
    GD.Weights{i,curlabel}     = Perf.Weights;
end
GD.DT{i,curlabel}          = Perf.TrDecisionValues;
GD.DV{i,curlabel}          = Perf.CVDecisionValues;
GD.C(i,:,curlabel)         = Perf.ModelComplexity.Complex;

if W2AVAIL
    GD.W2{i,curlabel}          = Perf.w2;
    GD.Md{i,curlabel}          = Perf.Md;
    GD.Mm{i,curlabel}          = Perf.Mm;
    GD.mMd(i,:,curlabel)       = mean(cellcat(GD.Md{i,curlabel}));
end

% Compute generalization error for current binary comparison
GD.ERR(i,:,curlabel)       = GD.TR(i,:,curlabel) - GD.TS(i,:,curlabel);

switch SVM.prog
    case 'SEQOPT'
        for curclass = 1:nclass
           GD.mSEQI{i, curclass, curlabel}      = Perf.MeanCritGain(curclass,:);
           GD.sdSEQI{i, curclass, curlabel}     = Perf.SDCritGain(curclass,:);
           GD.mSEQE{i, curclass, curlabel}      = Perf.MeanExamFreq(curclass,:);  
           GD.sdSEQE{i, curclass, curlabel}     = Perf.SDExamFreq(curclass,:); 
           GD.mSEQAbsThrU{i, curclass, curlabel} = Perf.MeanAbsThreshU(curclass,:);
           GD.sdSEQAbsThrU{i, curclass, curlabel} = Perf.SDAbsThreshU(curclass,:);
           GD.mSEQAbsThrL{i, curclass, curlabel} = Perf.MeanAbsThreshL(curclass,:);
           GD.sdSEQAbsThrL{i, curclass, curlabel} = Perf.SDAbsThreshL(curclass,:);
           GD.mSEQPercThrU{i, curclass, curlabel} = Perf.MeanPercThreshU(curclass,:);
           GD.sdSEQPercThrU{i, curclass, curlabel} = Perf.SDPercThreshU(curclass,:);
           GD.mSEQPercThrL{i, curclass, curlabel} = Perf.MeanPercThreshL(curclass,:);
           GD.sdSEQPercThrL{i, curclass, curlabel} = Perf.SDPercThreshL(curclass,:);
           GD.CasePropagations{i, curclass, curlabel} = nk_cellcat(CV2Perf.binCV1CasePropagations(:,:,curclass),[],2);
           GD.SeqPerfIncreases{i, curclass, curlabel} = nk_cellcat(CV2Perf.binCV1PerformanceIncreases(:,:,curclass),[],1);
           GD.DecValTraj{i,curclass,curlabel} = cat(3,CV2Perf.binCV1DecValTraj{:,:,curclass});
        end
    case 'WBLCOX'
        for curclass=1:nclass
           GD.mCutOffPerc(i, curclass, curlabel) = Perf.MeanThreshPerc(curclass,:);
           GD.sdCutOffPerc(i, curclass, curlabel) = Perf.SDThreshPerc(curclass,:);
           GD.mCutOffProb(i, curclass, curlabel) = Perf.MeanThreshProb(curclass,:);
           GD.sdCutOffProb(i, curclass, curlabel) = Perf.SDThreshProb(curclass,:);
           GD.CV1predictedtimes{i,curlabel} = Perf.PredictedTimes(:,:,curclass);
           GD.CV2predictedtimes{i,curlabel} = CV2Perf.binCV1times(:,:,curclass);
           GD.CV2Cutoffs{i,curlabel} = CV2Perf.binCV1probthresh(:,:,curclass);
        end            
end

%%%%% OPTIONALLY, PRINT INFO TO FIGURE:
% Create multi-group bar chart of current results

if ~BATCH && ~MULTILABEL.flag 

    if isfield(RAND,'CV2LCO') && ~isempty(RAND.CV2LCO), LCOstr = '_{LGO}'; end

    lefti = 0.10; wideelse = 0.15; sz = get( 0, 'Screensize' ); figuresz = sz;
    switch MODEFL
        case 'classification'
            for p = 1:DISP.nclass
                lg_long{p} = sprintf('M #%g: %s',p, CV.class{1,1}{p}.groupdesc);
                lg_short{p} = sprintf('M #%g',p);
            end
        case 'regression'
            lg_long{1} = sprintf('M #%g: %s',1,'Regression');
            lg_short{1} = sprintf('M #%g',1);
    end
    % Create position vectors depending on MULTI and flg
    if MULTI.flag
        figuresz(end) = 0.8*sz(end);
        height  = 0.35; bottom  = 0.47;
        bottomm = 0.05; 
    else
        figuresz(end) = 0.6*sz(end);
        height = 0.70; bottom = 0.1;
    end

    if flg
        figuresz(end-1) = 0.75*sz(end-1);
        widthperf = 0.4;
        pos1 = [lefti   bottom  widthperf    height];   % Performance axes
        pos2 = [0.60    bottom  wideelse     height];   % Complexity axes
        pos3 = [0.80    bottom  wideelse     height];   % Diversity axes
    else
        widthperf = 0.6;
        figuresz(end-1) = 0.5*sz(end-1);
        pos1 = [lefti   bottom  widthperf    height];   % Performance axes
        pos2 = [0.80    bottom  wideelse     height];   % Complexity axes
        pos3 = [];
    end

    if MULTI.flag
        posm1 = [ lefti   bottomm widthperf height ];  %
        if flg
            posm2 = [ 0.6 bottomm wideelse*2+0.05 height ];
        else
            posm2 = [];
        end
    end
    DISP.figuresz = figuresz;
    meanvec=zeros(1,nclass); stdvec=zeros(1,nclass); xaxlb=[]; xaxc=[];
    % Performance
    switch Metric
        case 1
            meanvec(1,:)     = Perf.MeanCVHTperf;
            stdvec(1,:)      = Perf.SDCVHTperf;
            xaxlb{1}         = 'Mean CV1 (T)';
        case 2
            meanvec(1,:)     = Perf.MeanCVHDperf;
            stdvec(1,:)      = Perf.SDCVHDperf;
            if SVM.RVMflag 
                xaxlb{end+1} = 'Mean CV1 (P)';
            else
                xaxlb{end+1} = 'Mean CV1 (D)';
            end
    end

    meanvec(end+1,:) = CV2Perf.binCV1Performance_Mean;
    stdvec(end+1,:)  = CV2Perf.binCV1Performance_SD;
    xaxlb{end+1}     = 'Mean CV2';

    % ****  Eventually include ensemble decisions *****
    switch Metric
        case 1
            meanvec(end+1,:) = CV2Perf.binCV2Performance_Targets;
            stdvec(end+1,:)  = 0;xaxlb{end+1} = 'Ensemble CV2 (T)';
        case 2
            meanvec(end+1,:) = CV2Perf.binCV2Performance_DecValues;
            stdvec(end+1,:)  = 0; 
            if SVM.RVMflag 
                xaxlb{end+1} = ['Ensemble CV2' LCOstr ' (P)'];
            else
                xaxlb{end+1} = ['Ensemble CV2' LCOstr ' (D)'];
            end
    end

    % Complexity
    meanc       = GD.C(i,:);
    stdc        = zeros(size(meanc));
    if DISP.nclass > 1
        meanc=[meanc;zeros(1,DISP.nclass)];
        stdc = [stdc;zeros(1,DISP.nclass)];
    end
    xaxc        = ' ';
    DISP.ax{1}.val_y = meanvec;
    DISP.ax{1}.std_y = stdvec;
    DISP.ax{1}.label = xaxlb;
    [DISP.ax{1}.ylm, DISP.ax{1}.ylb] = nk_GetScaleYAxisLabel(SVM);
    DISP.ax{1}.fc = 'b'; 
    DISP.ax{1}.title = 'Performance';
    DISP.ax{1}.position = pos1;
    DISP.ax{1}.lg = lg_long;

    DISP.ax{2}.val_y = meanc;
    DISP.ax{2}.std_y = stdc;
    DISP.ax{2}.label = xaxc;
    DISP.ax{2}.ylm = [0 100]; DISP.ax{2}.ylb = ' ';
    DISP.ax{2}.fc = 'r';
    DISP.ax{2}.title = 'Complexity';
    DISP.ax{2}.position = pos2;
    DISP.ax{2}.lg = [];

    if flg % Ensemble method
        meandiv = zeros(2,nclass); stddiv=zeros(2,nclass);
        meandiv(1,:)         = Perf.Ens_MeanCVDiv;
        meandiv(2,:)         = CV2Perf.binCV2Diversity_Targets;
        stddiv(1,:)          = Perf.Ens_SDCVDiv;
        xaxdiv             = {'Div CV1',         'Div CV2'};

        DISP.ax{3}.val_y = meandiv;
        DISP.ax{3}.std_y = stddiv;
        DISP.ax{3}.label = xaxdiv;
        DISP.ax{3}.ylm = [0 1]; 
        DISP.ax{3}.fc = 'g';
        if isfield(Param.EnsembleStrategy,'DivStr')
            DISP.ax{3}.ylb = Param.EnsembleStrategy.DivStr;
        else
            DISP.ax{3}.ylb = ' ';
        end
        DISP.ax{3}.title = 'Entropy';
        DISP.ax{3}.position = pos3;
        DISP.ax{3}.lg = [];
    end

    if MULTI.flag,
         %%%%% OPTIONALLY, PRINT INFO TO FIGURE:
         % Create multi-group bar chart of current results
         if ~MULTILABEL.flag

            DISP.m_ax{1}.val_y = [ Perf.MeanMultiTrPerf   Perf.MeanMultiCVPerf    mean(CV2Perf.MultiCV1Performance(:)) CV2Perf.MultiCV2Performance];
            DISP.m_ax{1}.std_y = [ Perf.SDMultiTrPerf     Perf.SDMultiCVPerf      std(CV2Perf.MultiCV1Performance(:))  0]; 
            DISP.m_ax{1}.label = {'Mean Tr',              'Mean CV1',             'Mean CV2',                          'Ensemble CV2'};
            DISP.m_ax{1}.ylm   = [0 100];
            if SVM.GridParam == 14
                DISP.m_ax{1}.ylb = 'Average Binary BAC [%]'; 
            else
                DISP.m_ax{1}.ylb = 'Mulit-group accuracy [%]'; 
            end
            DISP.m_ax{1}.title = 'Multi-group performance';
            DISP.m_ax{1}.position = posm1;
            DISP.m_ax{1}.lg = [];
            if flg
                DISP.m_ax{2}.val_y  = [GD.MultiM_DivT(i)         GD.MultiM_DivV(i)         GD.MultiCV2Div(i)];
                DISP.m_ax{2}.std_y  = [GD.MultiSD_DivT(i)        GD.MultiSD_DivV(i)        0]; 
                DISP.m_ax{2}.label  = {'Tr Div','CV1 Div', 'CV2 Div'};
                DISP.m_ax{2}.ylm    = [0 1];
                DISP.m_ax{2}.ylb    = 'Entropy'; 
                DISP.m_ax{2}.title = 'Multi-group diversity';
                DISP.m_ax{2}.position = posm2;
                DISP.m_ax{2}.lg = [];
            end
         end
     end
     DISP = nk_PrintCV2BarChart2(DISP, MULTI.flag);
end

if VERBOSE && ~MULTILABEL.flag
    for curclass=1:nclass% Print some info on the screen   
        fprintf('\nPerformance measures at parameter(s)')

        P_str = nk_DefineMLParamStr(DISP.P{curclass}, DISP.Pdesc, curclass);
        fprintf(' [%s]',P_str);

        switch MODEFL
            case 'classification'
                if RAND.Decompose ~=9
                     fprintf('\n%s: CV1 = %1.2f, CV2 = %1.2f, Complexity = %1.2f', CV.class{1,1}{curclass}.groupdesc, ...
                         GD.TR(i,curclass,curlabel), GD.TS(i,curclass,curlabel), GD.C(i,curclass,curlabel));
                else
                     fprintf('\nMulti-group: CV1 = %1.2f, CV2 = %1.2f, Complexity = %1.2f', ...
                     GD.TR(i,curclass,curlabel), GD.TS(i,curclass,curlabel), GD.C(i,curclass,curlabel));
                end
            case 'regression'
                 fprintf('\nRegression: CV1 = %1.2f, CV2 = %1.2f, Complexity = %1.2f', ...
                     GD.TR(i,curclass,curlabel), GD.TS(i,curclass,curlabel), GD.C(i,curclass,curlabel));
        end
    end

    if MULTI.flag
        fprintf('\nMulti-Group: CV1 = %1.2f, CV2 = %1.2f, Complexity = %1.2f', ...
            GD.MultiTR(i), GD.MultiTS(i), mean(GD.C(i,:)));
    end
end

% ___________________________________
function catarray = cellcat(arr, dim)

if nargin < 2, dim = 1; end
[ix,jx] = size(arr);
catarray = [];
for i=1:ix
    for j=1:jx
        switch dim 
            case 1
                catarray =[catarray; arr{i,j}];
            case 2
                catarray =[catarray arr{i,j}];
        end
    end
end

return