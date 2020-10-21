function [IO, act, mess] = DataIO(NM, parentstr, IO, mess, varind)
% =========================================================================
% function [IO, act, mess] = DataIO(NM, parentstr, IO, mess)
% -------------------------------------------------------------------------
% This function performs the main data import operations of NM
% It creates a flexible menu list of data input options and calls sanity
% check functions as well as the proper input functions for imaging data
% (ReadNifti, ReadSurf) or matrix-based data (ReadTabular). Reading
% text-based subject-level vector files has not been implemented yet.
%
% Input:
% ------
% NM            :       NeuroMiner structure
% parentstr     :       Name of calling parent function
% IO            :       Data input structure
% mess          :       Message structure containing sanity check results
% 
% Output:
% -------
% IO            :       Updated input structure
% act           :       Previously selected menu item
% mess          :       Updated message structure
% =========================================================================
% (c) Nikolaos Koutsouleris 06/2017

global SPMAVAIL FSAVAIL

oocvflag = 0;
modeflag = 'classification';
rowhead = 1;
colhead = 1;
n_samples = 0; 
na_str = '?';
spm_file = na_str;
spm_cov_avail = 0;
desc_groups{1} = na_str ; 
desc_data = na_str ;
brainmask = [];
badcoords = [];
case_edit = na_str ;
col_edit  = na_str ;
label_edit = na_str ;
matrix_edit = na_str;
globvar_edit = na_str;
labelmanage_edit = na_str;
covmanage_edit = na_str;
M_edit = na_str; 
gopt = na_str;
P = []; Pnan = [];
globscale = 1;
globnorm = 1;
delimit = 'tab'; sheet = na_str;
Thresh = struct('nVml',0,'Vml',0,'Lm',[],'threshop','gt');
completed = false;
reread_mat = false;
Yfl = 0;
wfu_flag = 2;
nangroup = 2; nan_subjects = 0;
if isfield(IO,'label') && sum(any(isnan(IO.label)))>0
    IO.nangroup = 1; IO.nan_subjects = sum(sum(isnan(IO.label),2));
end

labelmanage = {'add to existing','overwrite existing','skip'};
delimiters = {'comma','semi','space','tab','bar'};
globals = {'no globals','user-defined','computed by SPM'};
yn = {'yes','no'};
if ~exist('mess','var'), mess = []; end 

% Check if neuroimaging is available
if SPMAVAIL && FSAVAIL
    f = {'spm','nifti','surf','matrix'};
    datasource = 'nifti'; 
    groupmode = 0;
    
elseif SPMAVAIL
    f = {'spm','nifti','matrix'};
    datasource = 'nifti'; 
    groupmode = 0;
elseif FSAVAIL
    f = {'surf','matrix'};
    datasource = 'surf'; 
    groupmode = 0;
else
    f = {'matrix'};
    datasource = 'matrix';
    groupmode = 1;
end

% NM workspace structure containers
% Initialize IO
if ~exist('IO','var') || isempty(IO) || ~isstruct(IO)
    IO = struct('varind', varind, ...
               'datasource', datasource, ...
               'groupmode', groupmode, ...
               'spm_file', spm_file, ...
               'spm_cov_avail', spm_cov_avail, ...
               'labelmanage_edit', labelmanage_edit, ...
               'covmanage_edit', covmanage_edit, ...
               'nangroup', nangroup, ...
               'nan_subjects', nan_subjects, ...
               'M_edit', M_edit, ...
               'badcoords', badcoords, ...
               'brainmask', brainmask, ...
               'wfu_flag', wfu_flag, ...
               'globscale', globscale, ...
               'globnorm', globnorm, ...
               'globvar_edit', globvar_edit, ... 
               'matrix_edit', matrix_edit, ... 
               'label_edit', label_edit, ...
               'case_edit', case_edit, ...
               'col_edit', col_edit, ...
               'delimit', delimit, ...
               'sheet', sheet, ...
               'reread_mat', reread_mat, ...
               'completed', completed); 
end

if SPMAVAIL && (~isfield(IO,'spm_tbx') || isempty (IO.spm_tbx))
    IO.spm_tbx = spm('TBs');
    for i=1:numel(IO.spm_tbx)
        if strcmp(IO.spm_tbx(i).name,'wfupickatlas') 
            IO.wfu = true; 
            addpath(IO.spm_tbx(i).dir);
        end
        if strcmp(IO.spm_tbx(i).name,'cat12')
            IO.cat12 = true; 
            addpath(IO.spm_tbx(i).dir);
        end
    end 
end

% Check whether previous modality exists and copy respective pre-defined fields
% into IO to constrain the choices of the use
if isfield(NM,'Y') && ~isempty(NM.Y)
    Yfl = 1;
    if isfield(NM,'n_subjects')   
        IO.n_samples = numel(NM.n_subjects); 
        IO.n_subjects = NM.n_subjects;
        IO.n_subjects_all = NM.n_subjects_all;
    end
    if isfield(NM,'groupnames'),   IO.desc_groups = NM.groupnames; end
    if isfield(NM,'modeflag'),     IO.modeflag = NM.modeflag; end
    if isfield(NM,'label'),        IO.label = NM.label; end
    if isfield(NM,'covars'),       IO.covars = NM.covars; end
    if isfield(NM,'cases') && ~isfield(IO,'ID'), IO.ID = NM.cases; end
end

% Copy pre-existing settings of IO into active variables
if isfield(IO,'modeflag'),      modeflag = IO.modeflag; end
if isfield(IO,'n_samples'),     n_samples = IO.n_samples; end
if isfield(IO,'desc_groups'),   desc_groups = IO.desc_groups; end
if isfield(IO,'groupmode'),     groupmode = IO.groupmode; end
if isfield(IO,'desc_data'),     desc_data = IO.desc_data; end
if isfield(IO,'rowhead'),       rowhead = IO.rowhead; end
if isfield(IO,'colhead'),       colhead = IO.colhead; end
if isfield(IO,'brainmask'),     brainmask = IO.brainmask; end
if isfield(IO,'badcoords'),     badcoords = IO.badcoords;  end
if isfield(IO,'wfu_flag'),      wfu_flag = IO.wfu_flag;  end
if isfield(IO,'spm_cov_avail'), spm_cov_avail = IO.spm_cov_avail; end
if isfield(IO,'P'),             P = IO.P; end
if isfield(IO,'Pnan'),          Pnan = IO.Pnan; nan_subjects = IO.nan_subjects; end
if isfield(IO,'V'),             V = IO.V; end
if isfield(IO,'Thresh'),        Thresh = IO.Thresh; end
if isfield(IO,'globscale'),     globscale = IO.globscale;  end
if isfield(IO,'globnorm'),      globnorm = IO.globnorm;  end
if isfield(IO,'M_edit'),        M_edit = IO.M_edit; end
if isfield(IO,'matrix_edit'),   matrix_edit = IO.matrix_edit; end
if isfield(IO,'label_edit'),    label_edit = IO.label_edit;  end
if isfield(IO,'case_edit'),     case_edit = IO.case_edit;  end
if isfield(IO,'col_edit'),      col_edit = IO.col_edit; end
if isfield(IO,'labelmanage_edit'), labelmanage_edit = IO.labelmanage_edit; end
if isfield(IO,'covmanage_edit'), covmanage_edit = IO.covmanage_edit; end
if isfield(IO,'filt'),          filt = IO.filt; end
if isfield(IO,'filt_str'),      filt_str = IO.filt_str; end
if isfield(IO,'nangroup'),      nangroup = IO.nangroup; end
if isfield(IO,'modeflag'),      modeflag = IO.modeflag; else IO.modeflag = modeflag; end
if isfield(IO,'datasource'),    datasource = IO.datasource; end
if isfield(IO,'n_subjects'),    n_subjects = IO.n_subjects; end
if isfield(IO,'delimit'),       delimit = IO.delimit; end
if isfield(IO,'sheet'),         sheet = IO.sheet; end
if isfield(IO,'oocvflag'),      oocvflag = IO.oocvflag; else IO.labels_known=true; IO.oocvflag = false; end
if isfield(IO,'globvar_edit'),  globvar_edit = IO.globvar_edit ; end
if isfield(IO,'gopt'),          gopt = IO.gopt; else IO.gopt = gopt; end
if ~isfield(IO,'completed'),    IO.completed = completed; end
if ~isfield(IO,'reread_mat'),   IO.reread_mat = reread_mat; end

[IO, groupmode_str] = SetFileFilter(IO, groupmode, datasource);

mn_str = [];

%% Learning framework
if ~Yfl && ~oocvflag % NM workspace setup only
    mn_framework = sprintf('Define machine learning framework [ %s ]|',modeflag); mn_act = {'def_modeflag'};
else
    mn_framework = []; mn_act = {};
end
mn_str = [mn_str mn_framework]; 

%% Data provenience
if ~oocvflag, mn_filetype = sprintf('Select data input format [ %s ]|', datasource); mn_act = [ mn_act 'sel_dataorig' ]; mn_str = [mn_str mn_filetype]; end

%% Group definition (only for subject-level input files)
if ~Yfl && ~oocvflag % NM workspace setup only
    switch datasource 
        case 'spm'
            if isfield(IO,'design') && strcmp(IO.modeflag,'classification')
                groupstr = [];
                for i=1:numel(desc_groups)
                    groupstr = sprintf('%s, %s', groupstr, desc_groups{i}); 
                end
                mn_definegroups = sprintf('Modify sample identifiers [ %s ]|', groupstr(3:end));
                mn_act = [ mn_act 'def_groups' ]; mn_str = [mn_str mn_definegroups];
            end

        case {'nifti','surf','vector'}
            switch modeflag
                case 'classification'
                    if n_samples>0
                        groupstr = [];
                        for i=1:n_samples, groupstr = sprintf('%s, %s', groupstr, desc_groups{i}); end
                        mn_definegroups = sprintf('Modify no. samples [ N=%g ] and sample identifiers [ %s ]|', n_samples, groupstr(3:end));
                    else
                        mn_definegroups = 'Define no. of samples and sample identifiers [ ??? ]|';
                    end
                case 'regression'
                    regr_str = desc_groups{1};
                    mn_definegroups = sprintf('Provide sample identifier [ %s ]|', regr_str);                
            end
            mn_act = [ mn_act 'def_groups' ]; mn_str = [mn_str mn_definegroups];
            nangroup_str = yn{nangroup};
            mn_nangroup = sprintf('Involve unlabeled cases in the analyses [ %s ]|', nangroup_str);
            mn_act = [ mn_act 'def_nangroup' ]; mn_str = [ mn_str  mn_nangroup ];
            
        otherwise
            switch modeflag
                case 'regression'
                    regr_str = desc_groups{1};
                    mn_definegroups = sprintf('Provide sample identifier [ %s ]|', regr_str);
                    mn_act = [ mn_act 'def_groups' ]; mn_str = [mn_str mn_definegroups];
            end
    end
end

mn_groupmode_b = 'Define data provenance ';  matObj=[];
switch datasource
    
    %% Set neuroimaging-related input options
        
    case {'spm','nifti','surf'}
        
        if ~oocvflag
            
            if strcmp(datasource,'spm')
                [IO, mess] = SPMimport(IO,Thresh, 'check_matrix', mess);
                mn_spmpath = sprintf('Select path to SPM.mat [ %s ]|', IO.spm_file);
                mn_spmdesign = 'Review SPM design|';
                mn_act = [ mn_act 'sel_spmpath' ]; mn_str = [ mn_str mn_spmpath ];
                if ~strcmp(IO.spm_file,na_str)
                    mn_act = [ mn_act 'rev_spmdes' ]; mn_str = [ mn_str mn_spmdesign ];
                    if isfield(IO,'sel_dummy')
                        colvar_edit = 'Label Columns: '; colf = IO.sel_dummy;
                        for i=1:numel(colf), colvar_edit = sprintf('%s%g, ',colvar_edit, colf(i)); end
                        colvar_edit(end-1:end)=[];
                    else
                        colvar_edit = na_str;
                    end

                    % If previous label is found then the use can add new
                    % labels to the exsiting one, but only in the regression
                    % framework
                    if isfield(IO,'label') && strcmp(IO.modeflag,'regression')
                        if ~strcmp(IO.labelmanage_edit,na_str)
                            labelmanage_edit_str = labelmanage{IO.labelmanage_edit};
                        else
                            labelmanage_edit_str = na_str;
                        end
                        mn_managelabels = sprintf('Manage label extraction from SPM to NM workspace [ %s ]|', labelmanage_edit_str);
                        mn_act = [ mn_act 'sel_labelmanage' ]; mn_str = [ mn_str mn_managelabels ];
                    else
                        IO.labelmanage_edit = 2; % Do not add to existing labels
                    end

                    if ~isfield(IO,'label') || ~strcmp(IO.labelmanage_edit, na_str) && IO.labelmanage_edit < 3
                        mn_spmcolextract = sprintf('Specify which columns to extract as labels from SPM design [ %s ]|',colvar_edit);
                        mn_act = [ mn_act 'sel_spmlabelcols' ]; mn_str = [ mn_str mn_spmcolextract ];
                    end

                    if IO.spm_cov_avail

                         if isfield(IO,'covars')
                            if ~strcmp(IO.covmanage_edit,na_str)
                                covmanage_edit_str = labelmanage{IO.covmanage_edit};
                            else
                                covmanage_edit_str = na_str;
                            end
                            mn_managecov = sprintf('Manage covariate extraction from SPM to NM workspace [ %s ]|', covmanage_edit_str);
                            mn_act = [ mn_act 'sel_covmanage' ]; mn_str = [ mn_str mn_managecov ];
                         else
                            IO.covmanage_edit = 1; 
                        end
                        if ~isfield(IO,'covars') && ~strcmp(IO.covmanage_edit,na_str) && IO.covmanage_edit < 3
                            colvar_edit = 'Covariate Columns: '; 
                            if isfield(IO,'sel_dummy_covs')
                                if ~IO.sel_dummy_covs
                                    colvar_edit = sprintf('%s none',colvar_edit);
                                else
                                    colf = IO.sel_dummy_covs;
                                    for i=1:numel(colf), colvar_edit = sprintf('%s%g, ',colvar_edit, colf(i)); end
                                    colvar_edit(end-1:end)=[];
                                end
                            else
                                colvar_edit = na_str;
                            end
                            mn_spmcolextract = sprintf('Specify which columns to extract as covariates from SPM design [ %s ]|',colvar_edit);
                            mn_act = [ mn_act 'sel_spmcovcols' ]; mn_str = [ mn_str mn_spmcolextract ];
                        end
                    end
                end
            end
            
            if ~isempty(brainmask), [~,bmask_str] = fileparts(brainmask); else bmask_str = na_str; end
            if ~strcmp(bmask_str,na_str) && ~exist(brainmask,'file'), 
                bmask_act = 'Update';
                mess = GenerateMessageEntry(mess, sprintf('ERROR: Space-defining defining image could not be found in %s!\nUpdate path to point to the correct image!',brainmask));
            else
                bmask_act = 'Select';
            end
            if isempty(Thresh.Lm)
                mn_groupmode = sprintf('%s space-defining image [ %s ]|', bmask_act, bmask_str);
            else
                if numel(Thresh.Lm )> 1
                    mn_group_str = 'Parcellation atlas';
                elseif Thresh.nVml > 0
                    mn_group_str = 'Binarized image';
                else
                    mn_group_str = 'Binary mask';
                end
                mn_groupmode = sprintf('Select space-defining image [ %s: %s ]|', mn_group_str, brainmask);
            end
            mn_act = [mn_act 'def_spacedefimg']; mn_str = [mn_str mn_groupmode];
            if Thresh.nVml>0
                if numel(Thresh.Vml) == 1
                    mn_threshmode_str = sprintf('Threshold image at %g, Threshold operator: %s', Thresh.Vml, Thresh.threshop{1});
                else
                    mn_threshmode_str = sprintf('Extract %g unique values identified in space-defining image', numel(Thresh.Vml));
                end
                mn_threshmode = sprintf('Threshold space-defining image [ %s ]|', mn_threshmode_str);
                mn_act = [ mn_act 'sel_spacedefimgopt' ]; mn_str = [ mn_str mn_threshmode ];
            end
            if isfield(IO,'Pw') && ~isempty(IO.Pw)
                mn_weightimg_str = sprintf('Yes (%g images specified)',size(IO.Pw,1));
            else
                mn_weightimg_str = 'No';
            end
            mn_weightimg = sprintf('Read in additional image(s) for weighting [ %s ]|',mn_weightimg_str);
            mn_act = [ mn_act 'sel_weightimg' ]; mn_str = [ mn_str mn_weightimg ];
        else
            [~,bmask_str] = fileparts(brainmask); 
            if ~strcmp(bmask_str,na_str) && ~exist(brainmask,'file'), 
                bmask_act = 'Update';
                mess = GenerateMessageEntry(mess, sprintf('ERROR: Space-defining defining image could not be found in %s!\nUpdate path to point to the correct image!',brainmask));
                mn_groupmode = sprintf('%s space-defining image [ %s ]|', bmask_act, bmask_str);
                mn_act = [mn_act 'def_spacedefimg']; mn_str = [mn_str mn_groupmode];
            end
        end
        
        if n_samples && ~strcmp(brainmask, na_str) && ( strcmp(datasource,'nifti') || strcmp(datasource,'surf'))
            
            if ~isempty(P) && iscell(P) || (oocvflag && n_samples>1)
                filestr = []; mn_files = cell(n_samples,1);
                if ( ~oocvflag || ( oocvflag && IO.labels_known ))
                    for i=1:n_samples
                        if i>numel(P) || isempty(P{i})
                            filestri = 'NA';
                        else
                            filestri = num2str(size(P{i},1));
                        end
                        filestr = sprintf('%s, N(%s)=%s', filestr, desc_groups{i},filestri);
                        mn_files{i} = sprintf('Map image files to sample [ %s ]', filestr);   
                    end
                else
                    filestri = num2str(size(P{1},1));
                    filestr = sprintf('%s, N=%s',filestr, filestri);
                    mn_files = sprintf('Select sample [ %s ]', filestr);   
                    if (isfield(IO,'PP') && ~isempty(IO.PP)) && isfield(IO,'L')
                        nP=size(IO.PP,1); nL=size(IO.L,1);
                        %if IO.nangroup && IO.nan_subjects > 0; nL = nL + IO.nan_subjects; end
                        if nP~= nL, mess = GenerateMessageEntry(mess, sprintf('ERROR: No. of specified images (N=%g) does not match no. of cases (N=%g) in label array!',nP,nL)); end
                    end
                end
            else
                filestr = na_str;
                mn_files = sprintf('Select images files [ %s ]|', filestr);
            end
        
            if iscell(mn_files), 
                mn_files = strjoin(mn_files,'|');  mn_files_act = repmat('sel_img',n_samples,1) num2str(vec')]
            end
            mn_act = [ mn_act 'sel_img' ]; mn_str = [ mn_str mn_files ];
            
            if ~oocvflag && nangroup == 1
                if isempty(Pnan),
                    nanfilestr = na_str;
                else
                    nanfilestr = sprintf('N=%g',nan_subjects);
                end
                mn_nanfiles = sprintf('Select images files for the unlabeled cases [ %s ]|', nanfilestr);
                mn_act = [ mn_act 'sel_nanimg' ]; mn_str = [ mn_str mn_nanfiles ];
            end
            
            if isfield(IO,'Vinfo')
                mess = CheckImageInfoConsistency(IO.Vinfo, IO.Vvox, mess);
            end
            
            if strcmp(modeflag, 'regression') && (isfield(IO,'PP') && ~isempty(IO.PP)) && ~isfield(IO,'label')
                if ~oocvflag || IO.labels_known
                    Mx = 'MATLAB workspace'; groupmode_str = Mx;
                    groupmode_varstr = '(comma separated) variable name(s) (double)';
                    if ~strcmp(IO.label_edit,na_str), [ check_str, mess, IO.L ] = CheckLabel(IO, 'check_label', IO.PP, na_str, mess, Mx); else, check_str = na_str; end
                    mn_label_edit = sprintf('Define name of the MATLAB workspace variable containing the labels [ %s ]|', check_str); 
                    mn_act = [ mn_act 'def_labels' ]; mn_str = [mn_str mn_label_edit];
                    if isfield(IO,'L') && size(IO.L,2) > 1
                        if isfield(IO,'desc_labels')
                            desc_labels_str = strjoin(IO.desc_labels,', ');
                        else
                            desc_labels_str = na_str;
                        end
                        mn_labelnames_edit = sprintf('MULTI-LABEL MODE: Provide descriptions for these %g labels [ %s ]|', size(IO.L,2), desc_labels_str ); 
                        mn_act = [ mn_act 'def_label_names' ]; mn_str = [mn_str mn_labelnames_edit];
                    end
                end
            end
            
            if (isfield(IO,'PP') && ~isempty(IO.PP)) 
               
               if ~oocvflag, 
                   globscalestr = globals{globscale};
                   mn_globnorm = sprintf('Define global multiplier [ %g ]|',globnorm); mn_act = [ mn_act 'def_globnorm' ]; mn_str = [mn_str mn_globnorm]; 
                   mn_globscale = sprintf('Adjust data for globals [ %s ]|',globscalestr); mn_act = [ mn_act 'sel_globscale' ]; mn_str = [mn_str mn_globscale]; 
               end
               
               if globscale == 2
                    switch gopt
                        case 1
                            [ check_str, mess, IO ] = CheckGlobals(IO, 'check_globvar', na_str, mess); 
                            mn_globvar_edit = sprintf('Define name of the MATLAB workspace variable containing the globals vector [ %s ]|', check_str);
                            mess = GenerateMessageEntry(mess, 'The globals variable must consist of a vector containing the global values in the exactly the same order as the images.', 0, 'blue');
                        case 2
                            [ check_str, mess, IO ] = CheckGlobals(IO, 'check_globfile', na_str, mess);
                            mn_globvar_edit = sprintf('Provide path to globals file [ %s ]|', check_str); 
                            mess = GenerateMessageEntry(mess, sprintf(['The globals text file must contain two tab-delimited and header-less columns' ...
                                                       '\nwith rows ordered according to the order of the images:' ...
                                                       '\n\t1st column: Paths to %g image files', ...
                                                       '\n\t2nd column: Global volumes of %g cases'], size(IO.PP,1),size(IO.PP,1)), 0, 'blue');
                    end
                    mn_str = [ mn_str sprintf('%s',mn_globvar_edit) ]; mn_act = [ mn_act 'def_globvar' ];
               end
            end
        end
        if (isfield(IO,'PP') && ~isempty(IO.PP)) && strcmp(IO.datasource,'nifti')
            mn_display_images = 'Inspect image information and check registration|'; mn_act = [ mn_act 'disp_img' ]; mn_str = [mn_str mn_display_images ] ; 
        end
    %% Set matrix-related input options
    case 'matrix'
        flproc = false;
        switch modeflag
            case 'classification'
                mn_label_edit_m = 'cell array of strings';
                flproc = true;
            case 'regression'
                mn_label_edit_m = 'double vector/matrix';
                if ~strcmp(desc_groups{1},na_str), 
                    flproc=true; 
                else
                    mess = GenerateMessageEntry(mess, sprintf('Provide a sample identifier'));
                end        
        end
        
        if flproc
            
            switch groupmode 

                case {1,2} % MATLAB workspace or file as data source 

                    groupmode_varstr = 'variable name(s) (cell array of strings)'; 
                    switch groupmode
                        case 1 %%%%%% MATLAB workspace
                            groupmode_str = 'MATLAB workspace';

                        case 2 %%%%%% MATLAB file
                            groupmode_str = 'MATLAB mat file';
                            [p,f,e] = fileparts(IO.M_edit);
                            try 
                                matObj = matfile(IO.M_edit);
                            catch
                                mess = GenerateMessageEntry(mess, sprintf('ERROR: File ''%s%s'' not found in %s!',f,e,p));
                            end
                    end

                    mn_groupmode = sprintf('%s[ %s ]|',mn_groupmode_b, groupmode_str);
                    mn_act = [ mn_act 'sel_matsource' ]; mn_str = [mn_str mn_groupmode];

                    if ( groupmode == 2 && ~strcmp(M_edit, na_str) && ~strcmp(M_edit,'') ) || groupmode == 1

                        switch groupmode 
                            case 1
                                if ~strcmp(M_edit,na_str), [check_str, mess, IO] = CheckTabFile(IO, 'check_matrix', na_str, mess); else check_str = na_str; end
                                mn_matrix_edit = sprintf('Enter name of matrix variable containing the predictor data [ %s ]|', check_str);
                            case 2
                                if ~strcmp(matrix_edit,na_str),[check_str, mess, IO] = CheckTabFile(IO, 'check_matrix', na_str, mess); else check_str = na_str; end
                                mn_matrix_edit = sprintf('Enter name of matrix variable containing the predictor data [ %s ]|', check_str);
                        end
                        mn_act = [ mn_act 'def_matrix' ]; mn_str = [mn_str mn_matrix_edit];
                        if ~CheckReadyStatus(mess)
                            %% Collect info about and check the label array
                            if (~oocvflag || IO.labels_known) 
                                if ~strcmp(IO.label_edit,na_str), [ check_str, mess, IO ] = CheckTabFile(IO, 'check_label', na_str, mess); else check_str = na_str; end
                                mn_label_edit = sprintf('Enter name of label variable (%s) in %s [ %s ]|', mn_label_edit_m, groupmode_str, check_str); 
                                mn_act = [ mn_act 'def_labels' ]; mn_str = [mn_str mn_label_edit];
                            end

                            %% Collect info about and check the cases array
                            if ~strcmp(IO.case_edit,na_str), [ check_str, mess, IO] = CheckTabFile(IO, 'check_ID', na_str, mess); else check_str = na_str; end
                            mn_case_edit = sprintf('Enter name of case ID variable (cell array of strings ) in %s [ %s ]|', groupmode_str, check_str); 
                            mn_act = [ mn_act 'def_cases' ]; mn_str = [mn_str mn_case_edit];

                            %% Collect info about and check the features array
                            if ~strcmp(IO.col_edit,na_str), [ check_str, mess, IO ] = CheckTabFile(IO, 'check_featnames', na_str, mess); else check_str = na_str; end
                            mn_col_edit = sprintf('Enter name of feature descriptor variable (cell array of strings) in %s [ %s ]|', groupmode_str, check_str ); 
                            mn_act = [ mn_act 'def_featnames' ]; mn_str = [mn_str mn_col_edit];
                        end
                    end

                case {3,4} % Tabular text file or spreadsheet as source

                    groupmode_varstr = 'column header name(s)'; 
                    switch groupmode
                        case 3
                            groupmode_str = 'Text file'; 
                        case 4
                            groupmode_str = 'Spreadsheet file';
                    end

                    if isfield(IO,'M_edit') && ~strcmp(IO.M_edit,na_str) 

                        if ~strcmp(IO.M_edit,'')

                            switch groupmode 

                                case 3
                                    mn_groupmode = sprintf('%s [ Matrix from text file: %s ]|',mn_groupmode_b, M_edit);
                                    mn_act = [ mn_act 'sel_matsource' ]; 
                                    mn_str = [mn_str mn_groupmode];

                                    if ~oocvflag
                                        mn_delimiter = sprintf('Define delimiter in text file [ %s ]|',delimit); 
                                        mn_act = [ mn_act 'sel_delimit' ]; 
                                        mn_str = [ mn_str mn_delimiter ]; 
                                    end

                                case 4
                                    mn_groupmode = sprintf('%s [ Matrix from spreadsheet file: %s ]|',mn_groupmode_b, M_edit);
                                    mn_act = [ mn_act 'sel_matsource' ]; 
                                    mn_str = [mn_str mn_groupmode];
                                    if ~isfield(IO,'sheets')|| isempty(IO.sheets) || ~iscell(IO.sheets)
                                        fprintf('\nChecking spreadsheet file %s ...',IO.M_edit);
                                        [~,IO.sheets] = xlsfinfo(IO.M_edit);
                                    end
                                    if numel(IO.sheets)>1
                                        mn_sheet = sprintf('Define name of sheet in spreadsheet file [ %s ]|',sheet); mn_act = [ mn_act 'sel_sheet' ]; mn_str = [ mn_str mn_sheet ]; 
                                    else
                                        IO.sheet = IO.sheets{1};
                                    end
                            end

                            if ~oocvflag || IO.labels_known
                                [check_str, mess, IO] = CheckTabFile(IO, 'check_label', na_str, mess);
                                if isfield(IO,'t_Y') && ~isempty(IO.t_Y),cprintf('black*','Successfully opened %g x %g table file.',size(IO.t_Y,1),size(IO.t_Y,2)); end
                                mn_label_edit = sprintf('Specify column header containing the label data [ %s ]|', check_str); 
                                mn_act = [ mn_act 'def_labels' ]; mn_str = [mn_str mn_label_edit];
                            end

                            %% Collect info about and check the cases array
                            [check_str, mess, IO] = CheckTabFile(IO, 'check_ID', na_str, mess);
                            mn_case_edit = sprintf('Define column header containing the case IDs [ %s ]|', check_str); 
                            mn_act = [ mn_act 'def_cases' ]; mn_str = [mn_str mn_case_edit];

                            if isempty(mess) && ~strcmp(IO.case_edit,na_str) && ~any(strcmp(IO.label_edit,na_str)) 
                                [~, mess, IO] = CheckTabFile(IO, 'check_matrix', na_str, mess);
                            end
                        end
                    else
                        mn_groupmode = sprintf('%s[ %s: %s ]|',mn_groupmode_b, groupmode_str, na_str);
                        mn_act = [ mn_act 'sel_matsource' ]; mn_str = [mn_str mn_groupmode];
                    end

                case 5 %     XML read in (not implemented yet)

            end

            if (~CheckReadyStatus(mess) && (~sum(strcmp({M_edit,col_edit,case_edit},na_str)) || ~sum(strcmp({M_edit, matrix_edit, col_edit, case_edit},na_str))))
               mn_display_matrix = 'Inspect matrix data and select features for import|'; mn_act = [ mn_act 'm_inspect' ]; mn_str = [mn_str mn_display_matrix ] ; 
            end
        end

    case 'vector'
        
        % is structure (e.g. XML file)
        mn_vector_files = 'Select subject-level xml or csv files'; 
        if IO.vector.isstruct
            mn_vector_struct = sprintf('Define fieldname containing the data [ %s ]',IO.vector.structname); mn_act = [ mn_act 'def_structname' ]; mn_str = [ mn_str mn_vector_struct ];
        end
        mn_rowheader = sprintf('Define feature dimension [ %s ]|', IO.featuredim);  mn_act = [ mn_act 'def_rowhead' ]; mn_str = [ mn_str mn_rowheader ];
        mn_colheader = sprintf('Are variable names in the column headers [ %s ]|',mn_colheader_yesno);  mn_act = [ mn_act 'def_colhead' ];  mn_str = [ mn_str mn_colheader ];
        
end

if strcmp(modeflag,'regression') && ( any(~strcmp(IO.label_edit,na_str)) || (isfield(IO,'L') && ~isempty(IO.L))) 
    if ~oocvflag || IO.labels_known
        mn_showhist = 'Show histogram of label distribution|'; mn_act = [ mn_act 'disp_labelhisto' ]; mn_str = [ mn_str mn_showhist ];
    end
end

%% Data description
if ~oocvflag
    mn_descdata = sprintf('Describe data [ %s ]|', desc_data); 
    mn_act = [ mn_act 'def_desc' ]; mn_str = [mn_str mn_descdata];
end
% 
if strcmp(desc_data,na_str) && ~CheckReadyStatus(mess)
     mess = GenerateMessageEntry(mess,'* Describe the data modality to be generated!');
end

nk_PrintLogo; 

if ~isempty(matObj)
    fprintf('\n');cprintf('blue*','The %s contains the following variables: \n\n',groupmode_str); 
    whos(matObj)
end

disallow = CheckReadyStatus(mess);

if ~isempty(mess)
    for i=1:numel(mess)
        if isempty(mess(i).text), continue; end
        fprintf('\n');mess(i).text = regexprep(mess(i).text,'\','/');
        cprintf(mess(i).format,mess(i).text); 
    end
    fprintf('\n')
    mess = [];
end
if ~disallow
    mn_IO = sprintf('IMPORT %s',datasource ); mn_act = [ mn_act 'load_data' ] ; mn_str = [mn_str mn_IO];
end

fprintf('\n'); mestr = 'Input data into NM';  navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','You are here: %s >>> ',parentstr); 
act = char(nk_input(mestr,0,'mq', mn_str, mn_act));

switch act
    
    case 'BACK'
        IO.completed = false;
        return
    
    case 'def_modeflag'
        
        modeflag = nk_DefineModeflag(modeflag);
        if ~strcmp(IO.modeflag,modeflag)
            switch modeflag
                case 'regression'
                    % Collapse pre-specified samples into one sample when
                    % changing from classification to regression framework
                    if isfield(IO,'desc_groups') && numel(IO.desc_groups)>1
                        IO.desc_groups(2:end)=[]; 
                        if isfield(IO,'P')
                            IO.P{1} = IO.PP; IO.P(2:end)=[];
                            IO.V{1} = IO.VV; IO.V(2:end)=[];
                        end
                        IO.n_samples = 1;
                    end
                case 'classification'
                    if isfield(IO,'L'), IO = rmfield(IO,'L'); end
            end
        end
        IO.modeflag = modeflag;
        
    case 'sel_dataorig'
        % Select data provenance
        def = find(strcmp(f,datasource));
        if FSAVAIL
            fs_str = '|Cortical surface files (3D)';
        else
            fs_str='';
        end
        if SPMAVAIL
            spm_str = 'SPM (basic models)|ANALYZE/NIFTI images (3D)';
        else
            spm_str = '';
        end
        IO.datasource = char(nk_input('Select data origin',0,'m',sprintf('%s%s|Matrix data (2D)',spm_str,fs_str),f,def));
        switch IO.datasource
            case 'spm'
                IO.groupmode = -1;
                IO = SetFileFilter(IO, IO.groupmode, IO.datasource);
            case {'nifti','surf'}
                IO.groupmode = 0;
            otherwise
                IO.groupmode = 1;
        end
            
    case 'def_groups'
        switch datasource 
            case 'spm'
                IO.desc_groups = nk_input('Provide their names as a cell array of strings',0, 'e');
            case {'nifti','surf','vector','matrix'}
                [IO.n_samples, IO.n_subjects, IO.desc_groups, IO.labelflag] = nk_DefineGroups(oocvflag, modeflag);
        end
        
    case 'def_nangroup'
        IO.nangroup = nk_input('Involve unlabeled cases in the analysis',0,'yes|no',[1 2], nangroup);
        
    case 'sel_matsource'
        if ~oocvflag
            groupmode = nk_input('Select data matrix provenance',0,'m', ...
            'from MATLAB workspace|from MATLAB mat file|from a text file|from a spreadsheet file',1:4,groupmode);
        end
        IO = SetFileFilter(IO, groupmode, datasource);
        switch groupmode
            case {2,3,4,5}
                M_edit = nk_FileSelector(1,0,sprintf('Select %s file containing data, labels, case IDs, and feature names',groupmode_str), IO.filt);
        end
        switch groupmode
            case {3,4}
                IO.col_edit = ''; IO.matrix_edit = '';
        end
        IO.reread_mat = true;
        IO.groupmode = groupmode; IO.M_edit = M_edit;
        
    case 'sel_spmpath'
        IO.spm_file =  nk_FileSelector(1,0,'Select SPM.mat file', IO.filt);
        %[IO, mess] = SPMimport(IO,Thresh);
        
    case  'rev_spmdes'
        load(IO.spm_file)
        spm_DesRep('DesMtx',SPM.xX,SPM.xY.P); 
        spm_DesRep('DesRepUI',SPM);
        
    case 'sel_labelmanage'
        labelmanage_edit = nk_input('How should NM deal with possible labels in the SPM.mat',0,'mq', ...
            'Add to existing NM labels|Skip SPM label import',1:2,labelmanage_edit);
        if ~strcmp(labelmanage_edit,'BACK') || ~labelmanage_edit
            IO.labelmanage_edit = labelmanage_edit;
        end
        
    case 'def_label_names'
        if isfield(IO,'desc_labels')
            desc_labels = IO.desc_labels{i};
        else
            desc_labels = repmat({[]},size(IO.L,2),1);
        end
        for i=1:size(IO.L,2)
            IO.desc_labels{i} = nk_input(sprintf('Provide description for label #%g',i),0,'s',desc_labels{i});
        end
    case 'sel_covmanage'
        covmanage_edit = nk_input('How should NM deal with covariates in the SPM.mat',0,'mq', ...
            'Add to existing NM covariates|Overwrite existing NM covariates|Skip SPM covariate import',1:3,covmanage_edit);
        if ~strcmp(covmanage_edit,'BACK') || ~covmanage_edit
            IO.covmanage_edit = covmanage_edit;
        end
        
    case 'sel_spmlabelcols'
        if isfield(IO,'sel_dummy'), IO = rmfield(IO,'sel_dummy'); end
        [IO, mess ] = SPMimport(IO,Thresh,'check_label');
        
    case 'sel_spmcovcols'
        if isfield(IO,'sel_dummy_covs'), IO = rmfield(IO,'sel_dummy_covs'); end
        [IO, mess ] = SPMimport(IO,Thresh,'check_cov');
        
    case 'def_spacedefimg'
        if isfield(IO,'wfu') && IO.wfu && ~strcmp(IO.datasource,'surf') && ~oocvflag
            try
                wfu_str = sprintf('WFU Pickatlas (%s)',wfu_pickatlas_version);
            catch
                wfu_str = 'WFU Pickatlas';
            end
            IO.wfu_flag = nk_input('Retrieve space-defining image from',0,'mq',[wfu_str '|Space-defining image file'],[1,2],wfu_flag);
        else
            IO.wfu_flag = 2;
        end
        if IO.wfu_flag == 1
            IO.atlas_mask_filename = fullfile(pwd,sprintf('NM_mask_file_modality%g.nii',IO.varind));
            IO.wfu_atlas_region = wfu_pickatlas(IO.atlas_mask_filename);
            IO.brainmask = IO.atlas_mask_filename;
            IO.Vm = spm_vol(IO.brainmask);
            Lm = IO.wfu_atlas_region.names';
            IO.Thresh = nk_DataIO3_SpaceDefImage_config(IO.Vm, Lm);
            
        elseif IO.wfu_flag == 2
            [brainmask, Vm] = nk_FileSelector(1,IO.datasource,'Select space-defining image', IO.spacedef_filt);
            if isempty(brainmask), return, end;
            if isempty(Vm),
                mess = GenerateMessageEntry(mess,'ERROR: I can read MGZ files only in Linux!');
            elseif ~isempty(brainmask)
                if oocvflag
                    tVmO = IO.Vm; tVmO.fname = brainmask; tVmO = rmfield(tVmO,'private');
                    tVmN = Vm; tVmN = rmfield(tVmN,'private');
                    if ~isequal(tVmO, tVmN)
                        mess = GenerateMessageEntry(mess,'ERROR: The volume information of the updated space-defining image file does not match the old image information');
                    else
                        IO.brainmask = brainmask; IO.Vm = Vm; 
                    end
                else
                    IO.brainmask = brainmask; IO.Vm = Vm;
                    IO.badcoords = [];
                    IO.Thresh = nk_DataIO3_SpaceDefImage_config(IO.Vm, Thresh.Lm, IO.datasource);
                end
            end
        end
        
    case 'sel_spacedefimgopt'
        IO.Thresh = nk_DataIO3_SpaceDefImage_config(IO.Vm, Thresh.Lm);
        
    case 'sel_weightimg'
        if isfield(IO,'Pw') && ~isempty(IO.Pw), Pw=IO.Pw; useweight =1; else, Pw = []; useweight = 0; end 
        useweight = nk_input('Read-in images for weighting / subspace extraction',0,'yes|no',[1,0],useweight);
        if useweight,
            [Pw, Vw, mess] = nk_FileSelector(Inf, datasource, 'Select images for weighting / subspace extraction', IO.filt, Pw, mess);
        else
            Pw = []; Vw = [];
        end
        if ~isempty(Pw) && ~isempty(Vw), 
            IO.Pw = Pw; IO.Vw = Vw; 
        else
            if isfield(IO,'Pw'), IO = rmfield(IO,'Pw'); end
            if isfield(IO,'Vw'), IO = rmfield(IO,'Vw'); end
        end
        
    case 'sel_img'
        IO.PP = [];
        for i=1:n_samples,
            if ~oocvflag
                hdrstr = sprintf('Select %s for sample %g: %s', IO.datasource, i, desc_groups{i});
            else
                hdrstr = sprintf('Select %s for independent test data', IO.datasource );
            end
            if isfield(IO,'P') && numel(IO.P)>=i, Pi = IO.P{i}; else, Pi=[];  end
            [P, V, mess] = nk_FileSelector(n_subjects(i), datasource, hdrstr, IO.filt, Pi, mess);
            if ~isempty(P) && ~isempty(V)
                IO.P{i} = P; IO.V{i} = V; IO.PP = char(IO.PP,IO.P{i}); 
            else
                break
            end	
        end
        if ~isempty(IO.PP), 
            IO.PP(1,:)=[]; 
            [IO,mess] = RetrieveImageInfo(IO, datasource,mess);
            IO.n_subjects_all= size(char(IO.PP),1);
            [IO.ID, IO.files, mess] = nk_DefineCaseNames2(IO.PP, sum(IO.n_subjects_all),[], mess);
        end
        
    case 'sel_nanimg'
        hdrstr = sprintf('Select %s for unlabeled cases', IO.datasource );
        if nan_subjects > 0; ns = nan_subjects; else ns = Inf; end
        [IO.Pnan, IO.Vnan, mess] = nk_FileSelector(ns, datasource, hdrstr, IO.filt, mess);
        IO.nan_subjects = size(IO.Pnan,1);
        if ~isempty(IO.Pnan), 
            IO.PP = char(IO.PP,IO.Pnan);
            [IO,mess] = RetrieveImageInfo(IO, datasource,mess);
            IO.n_subjects_all= size(char(IO.PP),1);
            [IO.ID, IO.files, mess] = nk_DefineCaseNames2(IO.PP, sum(IO.n_subjects_all),[], mess);
        end
            
    case 'def_globnorm'
        IO.globnorm = nk_input('Specify global multiplier to be applied to all images',0,'e',globnorm);
    
    case 'sel_globscale'
        IO.globscale = nk_input('Proportionally scale all subject''s features to a global value (e.g. Total Intracranial Volume)',0,'m',...
                    'None|User-specified globals|Compute globals as mean voxel value',1:3,globscale);
        if ~isempty(P) && ~sum(cellfun(@isempty,P))
            switch IO.globscale
                case 1
                    IO.g = [];
                case 2  
                    IO.gopt = nk_input('Where can I find the global values',0,'mq','In MATLAB workspace variable|In text file',[1,2]);        
                case 3
                    IO = CalcGlobals(IO);
            end
        end
        
    case 'def_globvar'
        if ~isfield(IO,'n_subjects_all'), IO.n_subjects_all=sum(cellfun(@numel,IO.V));end
        switch gopt
            case 1
                IO.globvar_edit = nk_input(sprintf('Define name of globals variable [ '' %s '' ]', globvar_edit),0,'s',1);
            case 2
                IO.globvar_edit = nk_FileSelector(1,0,'Define path of globals text file', '.*\.txt$|.*\.dat$|.*\.csv');
        end
        
    case 'def_desc'
        IO.desc_data = nk_input(['Give a short description of the ' datasource ' data modality'],0,'s');
    
    case 'def_labels'
        if iscell(label_edit), label_edit = strjoin(label_edit,','); end
        label_edit = nk_input(sprintf('Define labels: %s in %s',groupmode_varstr, groupmode_str),0,'s', label_edit);
        t_label_edit = strsplit(label_edit,',');
        if numel(t_label_edit)==1,
            label_edit = deblank(char(t_label_edit));
        else
            label_edit = regexprep(t_label_edit,',','');
            label_edit = regexprep(label_edit,' ','');
        end
        IO.label_edit = label_edit;
        
    case 'def_cases'
        case_edit = nk_input(sprintf('Define cases: %s in %s', groupmode_varstr, groupmode_str),0,'s', case_edit);
       
        %if ~strcmp(case_edit,IO.case_edit), IO.reread_mat = true; end
        IO.case_edit = case_edit;
        
    case 'def_featnames'
        col_edit = nk_input(sprintf('Define name of the feature name variable (cell array of strings) in %s',groupmode_str),0,'s', col_edit);
        %if ~strcmp(col_edit,IO.col_edit), IO.reread_mat = true; end
        IO.col_edit = col_edit;
        
    case 'def_matrix'
        matrix_edit = nk_input(sprintf('Define name of data matrix in %s',groupmode_str),0,'s', matrix_edit);
        switch groupmode 
            case 1
                if ~strcmp(matrix_edit,IO.M_edit), IO.reread_mat = true; end
                IO.M_edit = matrix_edit;
            case 2
                if ~strcmp(matrix_edit,IO.matrix_edit), IO.reread_mat = true; end
                IO.matrix_edit = matrix_edit;
        end
        if ~oocvflag
            if isfield(IO,'selFeats'), IO = rmfield(IO,'selFeats'); end
        end
        if isfield(IO,'selCases'), IO = rmfield(IO,'selCases'); end
        
    case 'def_rowhead'
        IO.rowhead = nk_input('Are the case IDs in the row header?',0,'yes|no',[1,2], rowhead);
    
    case 'def_colhead'
        IO.colhead = nk_input('Are variable names in column headers?',0,'yes|no',[1,2], colhead);
    
    case 'sel_delimit'
        defdelimit = find(strcmp(delimiters,delimit));
        IO.delimit = char(nk_input('Specify column delimiter used in text file',0,'m','Comma|Semicolon|Space|Tabulator|Bar',{'comma','semi','space','tab','bar'},defdelimit)); 
        IO.reread_mat = true;
    
    case 'sel_sheet'
        mn_sel = []; mn_act = [];
        for i=1:numel(IO.sheets)
            mn_sel = sprintf('%s|%s',mn_sel,i,IO.sheets{i});
            mn_act = [mn_act i];
        end
        mn_sel(1)=[]; 
        sheet = IO.sheets{nk_input('Select the sheet you want to import the data from',0,'m',mn_sel,mn_act)};
        if ~strcmp(sheet,IO.sheet), IO.reread_mat = true; end
        IO.sheet = sheet;
        if isfield(IO,'selFeats'), IO = rmfield(IO,'selFeats'); end
        if isfield(IO,'selCases'), IO = rmfield(IO,'selCases'); end
        
    case 'disp_labelhisto'
        try
            t_label = evalin('base', IO.label_edit); 
        catch
            switch IO.groupmode
                case {-1,0,1}
                    t_label = IO.L;
                case 2
                    load(IO.M_edit, IO.label_edit); t_label = eval(IO.label_edit); 
                case 3
                    t_Y = readtable(IO.M_edit,'delimiter',IO.delimit);  
                    t_label = table2array(t_Y(:,{IO.label_edit}));
            end
        end
        figure;hist(t_label); 
        ylabel('# of measurements in bins'); 
        xlabel('Label measurements'); 
        title(sprintf('Histogram analysis of target label (in %s) for regression',groupmode_str)); 
        
 case 'disp_img'
        nk_PrintLogo
        if iscell(IO.PP), PP = char(IO.PP); F = char(IO.F); else, PP= IO.PP; F=IO.F; end
        fprintf('\n\n'); cprintf('black*','Found %g image files in setup:',size(PP,1));
        L = size(F,2); fnd = false(L,1);
        for i=1:size(F,1)
            mp = IO.Vinfo(i,:); mp_str = []; 
            if exist(deblank(PP(i,:)),'file')
                exist_str = 'File found'; fnd(i)=true;
                for j=1:3:12
                    j_str = sprintf('%g ',mp(j:j+2)); 
                    j_str(end)=[];
                    mp_str = [ mp_str j_str ' | '];
                end
                mp_str(end-2:end)=[]; 
            else
                exist_str = 'File NOT found!'; 
                mp_str = '';
            end
            cmdstr = ['\n%4g\t%' num2str(L) 's:\t%14s\t\t%s']; fprintf(cmdstr,i,F(i,:), exist_str, mp_str)
        end
        fprintf('\n');
        if any(fnd)
            img_ind = nk_input('Select images for display',0,'i',[1 size(PP,1)]);
            if isfield(IO,'Pw') && ~isempty(IO.Pw) 
                Pw=[]; for n=1:size(IO.Pw); if exist(deblank(IO.Pw(n,:)),'file'), Pw = char(Pw,IO.Pw(n,:)); end; end
                if ~isempty(Pw), 
                    Pw(1,:)=[]; PPsel = char(IO.brainmask, Pw, PP(img_ind,:)); 
                else
                    PPsel = char(IO.brainmask,PP(img_ind,:));
                end
            else
                PPsel = char(IO.brainmask,PP(img_ind,:));
            end
            spm_check_registration(PPsel);
        else
            nk_input('No image found. Press any key + enter to return',0,'sq'); 
        end
        
    case 'm_inspect'
        
        switch groupmode    
            case 1
                ID = evalin('base', case_edit); 
                F = evalin('base', col_edit);
                M = evalin('base',M_edit);
                
            case 2
                load(M_edit, case_edit); ID =eval(case_edit);
                load(M_edit, col_edit); F = eval(col_edit); 
                load(M_edit, matrix_edit); M =eval(IO.matrix_edit); 
                
            case {3,4}
                ID = table2cell(IO.t_Y(:,{case_edit}));
                F = IO.t_Y.Properties.VariableNames;
                col_cases = strcmp(F,case_edit);
                col_label = RetCellInCellsIndex(F, label_edit);
                M = table2array(IO.t_Y(:, ~col_cases & ~col_label ));
                F = F(~col_cases & ~col_label );
        end
        
        if isfield(IO,'selFeats'), selFeats = IO.selFeats; else, selFeats=[]; end
        if isfield(IO,'selCases'), selCases = IO.selCases; else, selCases=[]; end
        
        if isfield(IO,'label') 
            mode = 'feats'; 
        else
            if oocvflag
                mode = 'cases'; 
            else
                mode = 'all'; 
            end
        end
        
        [t_selFeats,t_selCases] = nk_ItemSelector('List', F, 'Matrix', M, 'Cases',ID, 'selFeats', selFeats, 'selCases', selCases, 'mode', mode);
        
        if ~isempty(t_selFeats), IO.selFeats = t_selFeats; end
        if ~isempty(t_selCases), IO.selCases = t_selCases; end
        
    case 'load_data'
        try   
            switch datasource 
                
                case {'spm','nifti'}
                    
                    IO.n_subjects = cellfun(@numel,IO.V); 
                    % Read NIFTI
                    IO = ReadNifti(IO);
                    if ~isfield(IO,'label')
                        % Read in label from workspace in case of regression
                        if strcmp(IO.modeflag,'regression') && strcmp(datasource,'nifti') && IO.labels_known
                            IO.L = evalin('base',IO.label_edit); 
                        end
                        %[IO.ID, IO.files] = nk_DefineCaseNames2(IO.PP, sum(IO.n_subjects_all));
                        IO = DefineLabels(IO, modeflag);
                    end
                    
                case 'surf'
                    
                    IO.n_subjects = cellfun(@numel,IO.V); 
                    % Read surface files
                    IO = ReadSurf(IO);
                    % Read in label from workspace in case of regression
                    if strcmp(IO.modeflag,'regression'), IO.L = evalin('base',IO.label_edit); end
                    IO = DefineLabels(IO, modeflag);
                   
                case 'matrix'
                    [IO, mess] = ReadTabular(IO, groupmode, mess);
                   
                case 'vector'
                    IO = ReadVectors(IO);
                    [IO.ID, IO.files] = nk_DefineCaseNames2(IO.PP, sum(IO.n_subjects_all));
                    IO = DefineLabels(IO);
            end
            flg = false;
            
            if ~isempty(mess),
                cnt = numel(mess)+1;
                for i=1:numel(mess)
                   if mess(i).flag == 1; flg = true; end
                end
            else
                cnt=1;
            end
            if flg 
                mess(cnt).format = 'red*'; 
                mess(cnt).text = sprintf(':-( Import of data didn''t work! Check your data import settings'); 
                mess(cnt).flag = 1;
                IO.completed = false;
            else
                mess(cnt).format = 'green*'; 
                mess(cnt).text = sprintf(':-D Import of data succeeded!'); 
                mess(cnt).flag = 0;
                IO.completed = true;
                act = 'BACK';
            end
        catch ERR
            mess.format = 'red*';
            mess.text = sprintf(['\n:-( Import of data didn''t work!!!' ...
                                 '\nMATLAB Error: %s (%s)' ...
                                 '\nFunction: %s (Line: %g)' ...
                                 '\nCheck your data import settings'], ...
                                 ERR.message, ERR.identifier,  ERR.stack(1).name, ERR.stack(1).line); 
            mess.flag = 1;
        end
        
        
end
