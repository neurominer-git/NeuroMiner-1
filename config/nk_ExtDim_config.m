function [dims, PercMode, act] = nk_ExtDim_config(RedMode, PercMode, dims, defaultsfl, parentstr)
global NM
if ~exist('defaultsfl','var') || isempty(defaultsfl); defaultsfl = false; end

if ~defaultsfl
    
    if (~exist('dims','var') || isempty(dims)) || (~exist('PercMode','var') || isempty(PercMode)); 
         [PercMode, dims] = nk_ExtDim_config(RedMode, [], [], 1);
    end
    defdims = dims; if size(defdims,2)==1 && numel(defdims)>1, defdims=defdims'; end
    mn_str = []; PercModeStr = {'Absolute range','Percentage range','Energy range'}; mn_act=[];
    switch RedMode
        case {'PCA', 't-SNE', 'SparsePCA'}
            mn_str = [mn_str sprintf('Define extaction mode for %s [ %s ]|',RedMode, PercModeStr{PercMode})]; mn_act = 1;
        case {'LDA', 'GDA'}
            PercMode = 1; L = NM.label; L(isnan(L))=[]; dims = numel(unique(L)); act=0; return
        otherwise
            PercMode = 1;
    end
    
    mn_str = [mn_str sprintf('Define extraction range [ %s ]',nk_ConcatParamstr(dims))]; mn_act = [mn_act 2];
    if numel(mn_act)>1
        nk_PrintLogo
        mestr = 'Extraction of components from reduced data projections'; navistr = [parentstr ' >>> ' mestr]; cprintf('*blue','\nYou are here: %s >>> ',parentstr); 
        act = nk_input(mestr,0,'mq',mn_str,mn_act);
    else
        act = 2;
    end
    switch act
        
        case 1
        
            PercMode = nk_input(sprintf('Define %s decomposition',RedMode),0,'m', ...
                            ['Absolute number range [ 1 ... n ] of eigenvectors|' ...
                             'Percentage range [ 0 ... 1 ] of max dimensionality|' ...
                             'Energy range [ 0 ... 1 ]of maximum decomposition'],1:3, PercMode);   
            switch PercMode 
                case {2,3}
                    dims = 0.8;
                case 1
                    dims = floor(size(NM.Y{NM.TrainParam.FUSION.M(1)},2)/10);
            end
        case 2
            switch PercMode
                case 1
                    inpstr = 'Dimensionalities to project data on (e.g: 1 5 10 or Start:Step:Stop)';
                case 2
                    inpstr = 'Percentages of max dimensionality as defined by training sample size (e.g: 25 50 75 or Start:Step:Stop)';
                case 3
                    inpstr = 'PCA ratios of complete decomposition (e.g: 0.25 0.5 0.75 or Start:Step:Stop)';
            end
            dims = nk_input(inpstr, 0, 'e', defdims); 
    end
     if numel(mn_act)<2, act = 0; end
else
    switch RedMode
        case {'PCA', 't-SNE'}
            dims = 0.8; PercMode = 3;
        case 'SparsePCA'
            dims = floor(numel(NM.cases)/10); PercMode = 1;
        case 'PLS'
            dims = 1;
        case {'LDA','GDA'}
            L = NM.label; L(isnan(L))=[]; dims = numel(unique(L)); 
        otherwise
            dims = floor(size(NM.Y{NM.TrainParam.FUSION.M(1)},2)/10); PercMode = 1;
    end
    act = 0;
end

end