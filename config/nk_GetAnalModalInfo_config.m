function inp = nk_GetAnalModalInfo_config(NM, inp)
global FUSION STACKING

if ~isfield(NM,'analysis'), inp=[]; return; end
if ~isfield(inp,'analind') || isempty(inp.analind)
    t_act = 1; while t_act>0, [t_act, inp.analind] = nk_SelectAnalysis(NM, 0, 'MAIN INTERFACE >> SELECT ANALYSIS', 1, [], 0, 0); end;
end
analysis = NM.analysis{inp.analind};

%% Define modality-independent parameters of current analysis as global variables
nk_SetupGlobVars2(analysis.params, 'setup_main', 0); 

if isfield(analysis,'rootdir')
    inp.rootdir = analysis.rootdir;
else
    inp.rootdir = analysis.GDdims{1}.RootPath;
end

if ~isempty(FUSION)        
    inp.F = analysis.params.TrainParam.FUSION.M; inp.nF = numel(inp.F); if FUSION.flag < 3, inp.nF = 1; end
else
    inp.F = 1; inp.nF = 1; FUSION.flag = 0; FUSION.M = 1;
end

inp.varind = FUSION.M; 
if ~isempty(STACKING) && STACKING.flag == 1
    inp.varind = 1;
    inp.varstr = '_varM';
    inp.metastr = '_META';
else
    inp.metastr = '';
    switch FUSION.flag
        case 0
            inp.varstr = sprintf('_var%g',FUSION.M); inp.concatfl=0;
        case 1
            inp.varind = 1;
            if numel(inp.F)>3
                inp.varstr = sprintf('_var_concat-%g',numel(inp.F));
            else
                inp.varstr = sprintf('_var%g',FUSION.M);
            end
        case 2
            inp.concatfl = 1; inp.varstr = [];
            %for i=1:numel(inp.varind), inp.varstr = sprintf('%s_var%g',inp.varstr,inp.varind(i)); end
        case 3
            inp.concatfl = 0; inp.varstr = [];
    end
end

