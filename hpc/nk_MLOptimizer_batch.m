function nk_MLOptimizer_batch(paramfile)

global NM

% Read parameter file
% ===========================
if ~exist(paramfile,'file') 
    error([paramfile 'not found. Abort job!']);
end
fid = fopen(paramfile);
params = textscan(fid, '%s');
fclose(fid);

% Parse parameters in params
% =======================================================================================================
NMpath        		= params{1}{1};					% NM root folder
datpath             = params{1}{2};					% NM structure
preprocmaster       = params{1}{3};                 % Optionally, path to PreprocMaster file
analind             = str2double(params{1}{4});		% Analysis index to identify analysis for HPC
ovrwrtfl            = str2double(params{1}{5}); 	% Runtime flag for MLOptimizer
curCPU              = str2double(params{1}{6});		% Current CPU to be used
numCPU              = str2double(params{1}{7});     % Number of CPUs in the submitted job
CV2x1               = str2double(params{1}{8});     % Range param for CV2 grid definition: Perm start CV2
CV2x2               = str2double(params{1}{9});     % Range param for CV2 grid definition: Perm end CV2
CV2y1               = str2double(params{1}{10});	% Range param for CV2 grid definition: Fold start CV2
CV2y2               = str2double(params{1}{11});	% Range param for CV2 grid definition: Fold end CV2

addpath(NMpath);
warning('off','MATLAB:FINITE:obsoleteFunction')
fprintf('\nLoading NM structure: %s',datpath)
load(datpath);
assignin('base','NM',NM);

% %%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE NeuroMiner %%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(preprocmaster) && exist(preprocmaster,'file')
    preprocmat = load(preprocmaster); lfl = 2; preprocmat = preprocmat.featmat;
else
    preprocmat = []; lfl = 1;
end

action = struct('addrootpath',1, ...
                'addDRpath',1, ...
                'addMIpath',1, ...
                'addLIBSVMpath',1, ...
                'addLIBLINpath',1, ...
                'addMikeRVMpath',1, ...
                'all',1);

nk_Initialize(action)

% %%%%%%%%%%%%%%%%%%%%%%%%% SETUP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
inp = struct('analind',			analind, ...
				'lfl',			lfl, ...
				'gdmat',		[], ...
				'gdanalmat', 	[], ...
				'varstr', 		[], ...
				'concatfl', 	[], ...
				'ovrwrt', 		ovrwrtfl, ...
				'update', 		true, ...
                'batchflag',    true);
				
inp.GridAct = nk_GenGridAct_batch(NM.analysis{analind}.params.cv, curCPU, numCPU, CV2x1, CV2x2, CV2y1, CV2y2);
inp.preprocmat = preprocmat;                            
inp = nk_GetAnalModalInfo_config(NM, inp);
                                
nk_MLOptimizerPrep(7, inp, 'NM:HPC:MLOPTIMIZER');
            
