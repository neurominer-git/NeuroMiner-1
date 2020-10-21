% =========================================================================
% FORMAT function [param, model] = nk_GetParam_LIBSVM(Y, label, ModelOnly, 
%                                                                 ...cmdstr)
% =========================================================================
% Train LIBSVM models and evaluate their performance using Y & label, 
% SlackParam, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2015

function [param, model] = nk_GetParam2_ROBSVM(Y, label, ModelOnly, Params)
                                            
global GRD EVALFUNC MODEFL CMDSTR VERBOSE                  

param = []; 
m = size(Y,1); if size(label,1) ~= size(Y,1), label = label'; end

% Check if sample weighting is necessary (currently regression only)
if strcmp(MODEFL,'regression') && SVM.LIBSVM.Weighting
    W = nk_WeigthDataInstanceHisto(label);
else
    W = ones(m,1);
end

GRDx = GRD;
GRDx.prog = 'ROBSVM'; CMDSTR.WeightFact=1;
% Check if weighting is necessary
cmd = nk_SetWeightStr(GRDx, MODEFL, CMDSTR, label, CMDSTR.cmd); %cmd = [cmd CMDSTR.quiet];

% Step 1: Get options and prepare for loops
options = nk_GenMatLearnOptions(Params);
options.corrmeth=0;
splitsz = ceil(m/options.nsplit);
%splitrem = m - options.nsplit*splitsz>0;
splitvec = 1:splitsz:m;
if splitvec(end)<m,splitvec(end+1)=m; end

D = zeros(m,1);

% Get test predictions
for i=1:numel(splitvec)-1
    trainind                                = true(m,1);
    trainind(splitvec(i):splitvec(i+1))     = false;
    testind                                 = ~trainind;
    model                                   = svmtrain312( W(trainind), label(trainind), Y(trainind,:), cmd );
    [~,~,D(testind)]                        = svmpredict312( label(testind), Y(testind,:), model );
end

% Identify errors
E = find(sign(D) ~= label);

% Standardize scores of error
Derr = D(E); Zerr = abs((Derr - median(Derr))./std(Derr));
remind = Zerr > options.wins;
trainind = true(m,1);
if sum(remind)>1
    Ex = E(remind); 
    switch options.corrmeth
        case 0 % Eliminate instance
            trainind(Ex) = false;
        case 1 % Flip label
            label(Ex) = sign(D(Ex));
    end
   if VERBOSE, fprintf('\t removed %g cases from training sample',sum(remind)); end
end

% Train full model with samples below threshold
model = svmtrain312( W(trainind), label(trainind), Y(trainind,:), cmd );

if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_LIBSVM([], Y, label, model) ;
    if GRD.ROBSVM.Optimization.b, param.dec_values = nk_CalibrateProbabilities(param.dec_values); end
    param.val = feval(EVALFUNC, label, param.target);
end

end