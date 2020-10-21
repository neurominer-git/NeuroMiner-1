function SUBGROUP = nk_SubgroupOperation_config( NM, SUBGROUP )

if isfield(NM,'C')
    menustr = ['Entire training set|' ...
               'User-defined subset of training partition|' ...
               'Random subset of the training population|' ...
               'Entire Calibration set|' ...
               'User-defined subset of the calibration set|' ...
               'Random subset of the calibation set']; menuact = 1:6;
else
    menustr = 'Entire training population|Subset of training population|Random subset of the training population'; menuact = 1:3;
end
SUBGROUP.flag = nk_input('Compute PCA modal from ',0,'m',menustr, menuact,SUBGROUP.flag); 
switch SUBGROUP.flag
    case 2
        SUBGROUP.ind = logical(nk_input('Define index vector (ones = used / zeros = unused) for subset computation',0,'e',[],[numel(NM.label),1]));
    case 3
        SUBGROUP.indperc = nk_input('Define percentage of population to be randomly selected',0,'i1',SUBGROUP.indperc) ;
    case 5
        SUBGROUP.ind = logical(nk_input('Define index vector (ones = used / zeros = unused) for subset computation',0,'e',[],[numel(NM.label),1]));
    case 6
        SUBGROUP.indperc = nk_input('Define percentage of population to be randomly selected',0,'i1',SUBGROUP.indperc) ;
end

