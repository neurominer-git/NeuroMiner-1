function dat = nk_ModelExportPrep(dat, analind)

global MODEFL CV FUSION MULTILABEL NM SAV

if isempty(NM), NM = evalin('base','NM'); end

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

dat.ExportParam.optfun  = 'CV1TS';
dat.ExportParam.optcrit = 'max';
dat.ExportParam.optperc = 95;
stranalysis             = SAV.matname;
inp1.PFE.Perc = 100;
inp1.PFE.Mode = 3;
% Loop through modalities
for i = 1:nF
    
    inp2 = nk_SetFusionMode2(dat, analysis, F, nF, i);
    inp = catstruct(inp1,inp2);
    
    for j = 1:MULTILABEL.dim
        
        if MULTILABEL.flag && MULTILABEL.dim>1
            fprintf('\n\n');cprintf('*black','====== Working on label #%g ====== ',j);
            inp.curlabel = j;
        else
            inp.curlabel = 1;
        end
        inp.multiflag = 0;
        inp.curmodal = i;
        
        strModelfile = fullfile(pwd,[stranalysis inp.varstr '_lb' num2str(j) '_ModelData_ID' dat.id '.mat']);
        
        try
            ModelData(i) = nk_ModelExport(inp);
            save(strModelfile,'ModelData');	
            dat.analysis{complvec(analind)}.ModelData.path = strModelfile;
            dat.analysis{complvec(analind)}.ModelData.status = 1;
        catch ERR
            dat.analysis{complvec(analind)}.ModelData.status = 0;
            errordlg(sprintf('NM was not able to export the model according to the chosen settings\n%s',ERR.message),'Model export problem');
        end
    end
end