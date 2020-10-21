function nk_ModelExport_config(dat, analind)

global SAV MODEFL CV OOCV FUSION MULTILABEL

if isempty(analind), return, end;
% Initialize analysis with current training parameters
complvec = [];
for z=1:numel(dat.analysis)
    if dat.analysis{z}.status, complvec = [ complvec z ]; end;
end
analysis = dat.analysis{complvec(analind)}; 
if ~isfield(analysis,'GDdims'), errordlg('No results detected for the selected analysis.\nRecreate analysis structure!'), end

% Define training parameters of current analysis as global variables
nk_SetupGlobVars2(analysis.params, 'setup_main', 0); 

if ~isempty(FUSION)        
    F = analysis.params.TrainParam.FUSION.M;
    nF = numel(F);
    if FUSION.flag < 3, nF = 1; end
else
    F = 1; nF = 1;
end

if strcmp(MODEFL,'classification')
    inp1.nclass = length(CV.class{1,1});
else
    inp1.nclass = 1;
end

% Loop through modalities
for i = 1:nF
    
    inp2 = nk_SetFusionMode2(dat, analysis, F, nF, i, oocvind);
    inp = catstruct(inp1,inp2);
    
    for j = 1:MULTILABEL.dim
        strModelfile = fullfile(pwd,[stranalysis inp.varstr '_lb' num2str(j) '_ModelData_ID' dat.id '.mat']);
        fprintf('\nLoading:\n%s',strOOCVfile);
		ModelData = nk_ModelExport(inp, GridAct);
        save(strModelfile,'ModelData');	
    end
end