function nk_Initialize(action)
% =========================================================================
% nk_Initialize(action)
% =========================================================================
% startup script for NeuroMiner 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2020

global NMinfo NM CALIBAVAIL OOCVAVAIL SPMAVAIL FSAVAIL

NMinfo.info.name    = 'NeuroMiner';
NMinfo.info.ver     = 'Release candidate version 1.05 | ELESSAR ';
NMinfo.info.author  = 'Nikolaos Koutsouleris';
NMinfo.info.affil   = 'Section for Neurodiagnostic Applications';
NMinfo.info.dep     = 'Department of Psychiatry and Psychotherapy';
NMinfo.info.inst    = 'Ludwig-Maximilian-University';
NMinfo.info.datever = '08/2020';
NMinfo.info.timestamp = date;
NMinfo.info.email   = 'nm@pronia.eu';
try
    NM = evalin('base','NM'); 
catch
    NM = [];
end
CALIBAVAIL          = false;
OOCVAVAIL           = false;

clc
fprintf('>>> Initializing NeuroMiner\n')
matpaths = path;

%if ~isdeployed
    if isunix, sep = ':'; else sep = ';'; end
    try
      matpaths = nk_strsplit(sep,matpaths);
    catch
      matpaths = nk_strsplit(matpaths, sep);
    end
    if exist(which('nm.m'),'file')
        neurominerpath = fileparts(which('nm.m'));
    else
        neurominerpath = fileparts(which('nm.p'));
    end
%else
    %neurominerpath = ctfroot;
%end

if nargin < 2
	action.all=1;
end

% Add root path
fprintf('\nRoot path is %s.\n',neurominerpath)
defs.path = neurominerpath; 
if action.addrootpath || action.all, addpath(defs.path); end

% Initialise SPM directory
SPMAVAIL = false; FSAVAIL=false;
imaging_init_path = fullfile(neurominerpath,'imaging_init.mat');
if ~isdeployed
    [spmrootdir, fsrootdir] = nk_ImagingInit(neurominerpath, imaging_init_path);
    if ischar(spmrootdir) && exist(spmrootdir,'dir') 
        SPMAVAIL = true;
        if ~checkpaths(matpaths,spmrootdir) 
            addpath(spmrootdir); 
            spm('Defaults','pet')
            fprintf('.');
        end  
    end
    if ischar(fsrootdir) && exist(fsrootdir,'dir')
        FSAVAIL = true;
        if ~checkpaths(matpaths,fsrootdir)
            addpath(fsrootdir); 
            fprintf('.');
        end  

    end
else
    if which('spm'), SPMAVAIL=true; spm('Defaults','pet'); FSAVAIL = true; end
end

% Initialize MEX file repository
if ~checkpaths(matpaths,fullfile(defs.path,'cfiles'))
    addpath(fullfile(defs.path,'cfiles'));
    fprintf('.');
end

if action.all
    
    %if ~isdeployed
        
        % Configuration subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'config'))
            addpath(fullfile(defs.path,'config'));
            fprintf('.');
            addpath(fullfile(defs.path,'config/data_import'));
            fprintf('.');
        end

        % Feature generation subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'FeatGen'))
            addpath(genpath(defs.path));
            addpath(fullfile(defs.path,'FeatGen/IMRelief')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/GFlip')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/Simba')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/RSS')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/RGS')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/FEAST')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/FEAST/FSToolbox')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/FEAST/MIToolbox')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/MRMR')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/SVM')); fprintf('.');
            addpath(fullfile(defs.path,'FeatGen/SA')); fprintf('.');
        end
        
        % Preprocessing subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'preproc'))
            addpath(fullfile(defs.path,'preproc')); %fprintf('.');
            addpath(fullfile(defs.path,'preproc/mi')); fprintf('.');
            addpath(fullfile(defs.path,'preproc/combat')); fprintf('.');
            addpath(fullfile(defs.path,'preproc/drtoolbox')); 
            addpath(fullfile(defs.path,['preproc/drtoolbox' filesep 'techniques'])); fprintf('.');
            addpath(fullfile(defs.path,'preproc/drtools'));fprintf('.');
            addpath(fullfile(defs.path,'preproc/drrobust'));fprintf('.');
            addpath(fullfile(defs.path,'preproc/nmfv1_4')); fprintf('.');
            addpath(fullfile(defs.path,'preproc/NeNMF')); fprintf('.');
            addpath(fullfile(defs.path,'preproc/NMF_Soteiras')); fprintf('.');
            addpath(fullfile(defs.path,'preproc/SparsePLS')); fprintf('.');
        end
        
        % Train & Predict functions subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'trainpredict'))
            addpath(fullfile(defs.path,'trainpredict')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/libsvm-weights-3.12/matlab')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/libsvm-mat-2.9-1')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/liblinear-2.20/')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/LIBSVM-Plus-2.89')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/mRVMs')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/MRVR'));fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/SB2_Release_200'));fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/elm_linear'));fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/blogreg')); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/matLearn2016'));
            addpath(genpath(fullfile(defs.path,'trainpredict/matLearn2016'))); fprintf('.');
            addpath(genpath(fullfile(defs.path,'trainpredict/glmnet'))); fprintf('.');
            addpath(fullfile(defs.path,'trainpredict/weibull_cox_v_1_1_4')); fprintf('.');
        end
        
        if ~checkpaths(matpaths,fullfile(defs.path,'gridsearch'))
            addpath(fullfile(defs.path,'gridsearch')); fprintf('.');
        end
        
        % Visualization subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'visual'))
            addpath(fullfile(defs.path,'visual')); fprintf('.');
        end
        
        % Visualization subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'gui'))
            addpath(fullfile(defs.path,'gui')); fprintf('.');
            addpath(fullfile(defs.path,'gui/export_fig')); fprintf('.');
        end

        % Ensemble learning subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'ensemble'))
            addpath(fullfile(defs.path,'ensemble'));fprintf('.');
        end

         % Utilities subdirectory
        if ~checkpaths(matpaths,fullfile(defs.path,'util'))
            addpath(fullfile(defs.path,'util')); fprintf('.');
            addpath(fullfile(defs.path,'util/freesurfer')); fprintf('.');
        end
    %end
end

if action.addDRpath || action.all
    
    %if ~isdeployed
        % Preprocessing subdirectory
        
        % Path of dimensionality reduction toolbox
        if ~checkpaths(matpaths,fullfile(defs.path,'preproc/drtoolbox'))
            addpath(fullfile(defs.path,'preproc/drtoolbox')); 
            addpath(fullfile(defs.path,['preproc/drtoolbox' filesep 'techniques'])); fprintf('.');     
        end

        if ~checkpaths(matpaths,fullfile(defs.path,'preproc/drtools'))
            addpath(fullfile(defs.path,'preproc/drtools')); fprintf('.');
        end

        if ~checkpaths(matpaths,fullfile(defs.path,'preproc/drrobust'))
            addpath(fullfile(defs.path,'preproc/drrobust')); fprintf('.');
        end
        
        % Path of dimensionality reduction toolbox
        if ~checkpaths(matpaths,fullfile(defs.path,'preproc/nmfv1_4'))
            addpath(fullfile(defs.path,'preproc/nmfv1_4')); fprintf('.');
        end
    
    %end
	
end
        
if action.addMIpath
    %if ~isdeployed
        if ~checkpaths(matpaths,fullfile(defs.path,'preproc/mi'))
            addpath(fullfile(defs.path,'preproc/mi')); fprintf('.');
        end
    %end
end

if action.addLIBSVMpath
    
    %if ~isdeployed
        if ~checkpaths(matpaths,fullfile(defs.path,'trainpredict/libsvm-weights-3.12/matlab'))
            addpath(fullfile(defs.path,'trainpredict/libsvm-weights-3.12/matlab')); fprintf('.');
        end
    %end
    
    %if ~isdeployed
        if ~checkpaths(matpaths,fullfile(defs.path,'trainpredict/libsvm-mat-2.9-1'))
            addpath(fullfile(defs.path,'trainpredict/libsvm-mat-2.9-1'));fprintf('.');
        end
    %end

    %if ~isdeployed
        if ~checkpaths(matpaths,fullfile(defs.path,'trainpredict/LIBSVM-Plus-2.89'))
            addpath(fullfile(defs.path,'trainpredict/LIBSVM-Plus-2.89'));fprintf('.');
        end
    %end

end
if action.addLIBLINpath 
    %if ~isdeployed 
        if ~checkpaths(matpaths,fullfile(defs.path,'trainpredict/liblinear-2.1/matlab'))
           addpath(fullfile(defs.path,'trainpredict/liblinear-2.1/matlab')); fprintf('.');
        end
    %end
end

if action.addMikeRVMpath 
    %if ~isdeployed 
        if ~checkpaths(matpaths,fullfile(defs.path,'trainpredict/SB2_Release_200'))
            addpath(fullfile(defs.path,'trainpredict/SB2_Release_200'));fprintf('.');
        end
    %end
end

if ~isempty(NM) && ~isstruct(NM) 
    error('The NM workspace variable does not have the correct format and thus cannot be not recognized. Clear the variable and start again!')
end
NM.defs.paths = defs.path;
NM.defs.NM_ver = NMinfo.info.ver;
if ~isfield(NM.defs,'ver'), NM.defs.ver = ver; end
if ~isfield(NM.defs,'computer'), NM.defs.computer = computer; end

% _______________________________________
function [fnd, pathx] = checkpaths(matpaths,pathx,findmode)

fnd = 0;
if ~exist('findmode', 'var') || isempty(findmode)
    findmode = 1;
end
for i=1:length(matpaths)
    switch findmode
        case 1
            if strcmp(matpaths{i}, pathx)
                pathx = matpaths{i};
                fnd=1;
                break
            end
        case 2
            if strfind(matpaths{1}, pathx)
                pathx = matpaths{i};
                fnd=1;
                break
            end
    end
end

return