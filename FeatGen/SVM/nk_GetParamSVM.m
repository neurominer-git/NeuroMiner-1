% =========================================================================
% FORMAT function [param, model] = nk_GetParamSVM(Y, label, SlackParam, ...
%                                                KernParam, ModelOnly)
% =========================================================================
% Train LIBSVM models and evaluate their performance using Y & label, 
% SlackParam, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 02/2012

function [param, model] = nk_GetParamSVM(Y, label, cmdstr, ...
                                            SVM, LIBSVMTRAIN, LIBSVMPREDICT, EVALFUNC, ...
                                            TrainParams, ModelOnly)
                                        
global MODEFL
                                                                
param = []; flw = 0;
if strcmp(LIBSVMTRAIN,'svmtrain312'), flw = 1; end

if isfield(TrainParams,'SlackParam')
    cmdstr = [cmdstr.simplemodel ' -c ' TrainParams.SlackParam ];
else
    cmdstr = [cmdstr.simplemodel ' -n ' TrainParams.NuCParam ];
end

CMDSTR.WeightFact = 1;
CMDSTR.quiet = 1;

% Check if weighting is necessary
cmdstr = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmdstr); 

if isfield(TrainParams,'NuParam'), cmdstr = [cmdstr ' -n ' TrainParams.NuParam]; end
if isfield(TrainParams,'EpsParam'), cmdstr = [cmdstr ' -p ' TrainParams.EpsParam]; end

cmdstr = [cmdstr CMDSTR.quiet];

% Check if sample weighting is necessary (currently regression only)
if flw
    W = ones(numel(label),1);
    if strcmp(MODEFL,'regression') && SVM.LIBSVM.Weighting
        W = nk_WeigthDataInstanceHisto(label);
    end
end

if iscell(Y) 
   
    % MKL-based learning not implemented yet
   
else % Univariate case
    if flw
         model = feval( LIBSVMTRAIN, W, label, Y, cmdstr );
    else
        model = feval( LIBSVMTRAIN, label, Y, cmdstr );
    end
    if ~ModelOnly
        [param.target, param.dec_values] = nk_GetTestPerfSVM(Y, label, model, SVM, LIBSVMPREDICT) ;
        param.val = feval(EVALFUNC, label, param.target);
    end
end