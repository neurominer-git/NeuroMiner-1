function param = nk_CostFun_config(param, SVM, MODEFL, defaultsfl)

% Defaults
if strcmp(SVM.kernel.kernstr,' -s 0') 
    CostFun = 3;
else
    CostFun = 2;
end

if isfield(param,'CostFun'), CostFun = param.CostFun; end
if ~exist('defaultsfl','var'), defaultsfl = 0; end;

if ~defaultsfl
    
   CostFun= nk_input('Cost function',0,'m', ...
    ['Optimize ensemble using CV1 training data ' MODEFL ' performance|' ...
    'Optimize ensemble using CV1 test data ' MODEFL ' performance|', ...
    'Optimize ensemble using CV1 training & CV1 test data ' MODEFL ' performance'],1:3,CostFun);
   
end

param.CostFun = CostFun;