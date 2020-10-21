function param = nk_SubSpaceStrategy_config(param, SVM, MODEFL, defaultfl, parentstr)
%
% Defaults
% --------
SubSpaceStrategy = 1;   % Winner takes it all
SubSpaceCrit = 0;       % No subspace threshold

% Take old values if available
if isfield(param,'SubSpaceStrategy'), SubSpaceStrategy = param.SubSpaceStrategy; end
if isfield(param,'SubSpaceCrit'), SubSpaceCrit = param.SubSpaceCrit; end

if ~exist('defaultfl','var') || isempty(defaultfl),  defaultfl = 0; end;

if ~defaultfl
    
    if ~param.CostFun, param.CostFun = 2; end
    d = nk_GetParamDescription2([],param,'EnsType');
    mestr = 'Subspace selection setup'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>>',parentstr); 
    
    nk_PrintLogo
    switch SVM.GridParam
        case {9,11,12,18}
            df = 25; mxstr ='minimum ';
            SubSpaceStrategy = nk_input(mestr,0, ...
            'm',['Subspace with ' d.CostType ' criterion (winner takes it all)|' ...
                'Subspaces with ' d.CostType ' within a range of the maximum|' ...
                'Subspaces with ' d.CostType ' BELOW a ? percentile|' ...
                'All Subspaces'],1:4, SubSpaceStrategy);

        otherwise
            df = 75; mxstr='maximum';
            SubSpaceStrategy = nk_input(mestr,0, ...
            'm',['Subspace with ' d.CostType ' criterion (winner takes it all)|' ...
                'Subspaces with ' d.CostType ' within a range of the maximum|' ...
                'Subspaces with ' d.CostType ' ABOVE a ? percentile|' ...
                'All Subspaces'],1:4, SubSpaceStrategy);

    end

    switch SubSpaceStrategy
        case {1,4}
            SubSpaceCrit = 0;
        case 2
            SubSpaceCrit = nk_input(['Range [%] from ' d.CostType mxstr],0,'e', 5);
        case 3
            SubSpaceCrit = nk_input(['Percentile [%] for ' d.CostType ' cutoff'],0, 'e', df);
    end
        
   %param = nk_EnsembleStrategy2_config(param, SVM, MODEFL, [], navistr);
    
else
   %param = nk_EnsembleStrategy2_config(param, SVM, MODEFL, 1);
end

param.SubSpaceStrategy              = SubSpaceStrategy;
param.SubSpaceCrit                  = SubSpaceCrit;
