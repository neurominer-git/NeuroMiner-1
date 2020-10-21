% =========================================================================
% FORMAT param = AUC(expected, predicted)
% =========================================================================
% Compute Area-Under-the-Curve for binary classification problems
% Based on the code of Chih-Jen Lin
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2015

function param = AUC(expected, predicted)
param = []; 
if isempty(expected), return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);
uE = unique(expected);
if numel(uE)==2, 
      param = fastAUC(expected, predicted,1);
%     [~,idx]         = sort(predicted, 'descend');
%     expected        = expected(idx);
%     tp              = cumsum(expected == 1);
%     fp              = sum(expected  == -1);
%     ret             = sum(tp(expected == -1));
%     if tp == 0 | fp == 0;
%         warning('Too few postive true labels or negative true labels');
%         param = 0;
%     else
%         param = ret / tp(end) / fp;
%     end
else
    error('AUC works only for binary problems'); 
end