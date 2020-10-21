% =========================================================================
% FORMAT function [param, model] = nk_GetParam_MikRVM(Y, label, ~, KernParam, 
%                                                                  ModelOnly)
% =========================================================================
% Train MikRVM models and evaluate their performance using Y & label, KernParam
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06/2011

function [param, model] = nk_GetParam2_MKLRVM(Y, label, ModelOnly, KernParam)

global SVM MKLRVM EVALFUNC

% Generate Label Matrix
if iscell(label)
    tlabel = label{1};
else
    tlabel = label;
end

t = nk_MakeRVMTargetMatrix(tlabel)';

if iscell(Y), 
    funcname = [MKLRVM.funcname_learn '_MKL']; 
else
    funcname = MKLRVM.funcname_learn;
end;

model = feval(funcname, '-p', Y, t, MKLRVM.standardize_flag, ...
    MKLRVM.converg, MKLRVM.nmax, SVM.kernel.kernstr, KernParam, ...
    MKLRVM.plot_flag);

if ~ModelOnly
    if iscell(Y), 
        funceval = [MKLRVM.funcname_predict '_MKL'];
    else
        funceval = MKLRVM.funcname_predict; 
    end;
    % Get sample group membership probabilities
    param.dec_values = feval(funceval, model, Y, t, Y);
    % Convert to predicted group memberships (binary case)
    [mx,param.target] = max(param.dec_values, [], 2); param.target(param.target == 2) = -1;
    param.dec_values = param.dec_values(:,1);
    % Compute performance of MKLRVM
    param.val = feval(EVALFUNC, label, param.target);
else
    param = [];
end
        