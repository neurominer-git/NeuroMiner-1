% =========================================================================
% FORMAT function [param, model] = nk_GetParam_kNNMEX(Y, label, ModelOnly, Params)
% =========================================================================
% Train kNN model using specified number (Params) of nearest neighbors 
% if ModelOnly = 1, return only model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2015

function [param, model] = nk_GetParam2_kNNMEX(Y, label, ModelOnly, Params)
global EVALFUNC

uL = unique(label); if uL(1)<0, model.sfl = true; label(label<0)=2; else, model.sfl = false; end
model.X = Y; model.Y = label; model.kNN = Params;
if ~ModelOnly
    [target, dec_values] = fknn(model.X, model.Y, model.X, [], model.kNN, 0);
    if model.sfl; target(target==2)=-1; end
    param.target = target; param.dec_values = dec_values(:,1);
    param.val = feval(EVALFUNC, label, param.target);
else
    param = [];
end