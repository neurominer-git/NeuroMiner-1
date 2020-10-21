% =========================================================================
% FORMAT param = nk_OOCV_config(param, res)
% =========================================================================
% NeuroMiner menu configurator for the OOCV runtime parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08 / 2020

function [param, act] = nk_OOCV_config(param)
global NM

meanflag = 1;
savemodels = 2;
saveoocvdata = 1;
if isfield(param,'MULTI') && param.MULTI.flag, 
    multiflag = true; groupmode = 2; 
else
    multiflag = false;
    groupmode = 1; 
end
trainwithCV2Ts = 1;

if ~isfield(NM.TrainParam,'OOCV')
    NM.TrainParam.OOCV.meanflag = meanflag;
    NM.TrainParam.OOCV.groupmode = groupmode;
    NM.TrainParam.OOCV.multiflag = multiflag;
    NM.TrainParam.OOCV.trainwithCV2Ts = trainwithCV2Ts;
end
o = nk_GetParamDescription2(NM,NM.TrainParam,'oocv');

nk_PrintLogo

if multiflag && ~strcmp(NM,'regression')
    act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> OOCV settings',0,'mq', ...
                    ['Ensemble mode [ ' o.meanflagstr ' ]|' ...
                     'Group processing mode [ ' o.groupmodestr ' ]|' ...
                     'Model retraining mode [ ' o.trainwithCV2Tsstr ' ]'],1:3);
else
    act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> OOCV settings',0,'mq', ...
                    ['Ensemble mode [ ' o.meanflagstr ' ]|' ...
                     'Model retraining mode [ ' o.trainwithCV2Tsstr ' ]'],[1 3]);
end

switch act

    case 1
        meanflag = nk_input('Ensemble construction mode', 0, 'm', ...
                                ['Aggregate ALL base learners into ONE ensemble|' ...
                                 'Compute mean decisions of CV1 partion ensembles before aggregating'], ...
                                 [1,2], meanflag);
    case 2
        groupmode = nk_input('Define group processing mode',0, 'm', ...
                                ['OOCV prediction at binary predictors'' optimum parameters|' ...
                                 'OOCV prediction at multi-group predictors'' optimum parameters|' ...
                                 'OOCV prediction at binary & multi-group predictors'' optimum parameters'], ...
                                 1:3, groupmode);
            
    case 3
        trainwithCV2Ts = nk_input('Model retraining mode', 0, 'm', ...
                                 ['As defined for CV training phase|' ...
                                  'Use all available CV data (CV1+CV2)'], ...
                                  [1,2], trainwithCV2Ts); 
end

param.OOCV.meanflag         = meanflag;
param.OOCV.groupmode        = groupmode;
param.OOCV.trainwithCV2Ts   = trainwithCV2Ts;
param.OOCV.savemodels       = savemodels;
param.OOCV.saveoocvdata     = saveoocvdata;

end