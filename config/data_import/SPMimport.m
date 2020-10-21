function [IO, mess] = SPMimport(IO, Thresh, act, mess)

if ~exist('mess','var'), mess=[]; end
if ~exist('ShowDesign','var') || isempty(ShowDesign), ShowDesign = true; end
    
try
    
    if ~exist(IO.spm_file,'file'), return, end
    
    load(IO.spm_file)
    
    switch act
            
        case 'check_matrix'
            
            %% Read files
            if ~isfield(IO,'PP') || isempty(IO.PP)
                IO.PP = regexprep(SPM.xY.P,',1','');
            end
            %% Read core SPM design elements
            IO.xX = SPM.xX; % Design matrix
            IO.xC = SPM.xC; % Covariates
            if ~isempty(IO.xC)
                IO.spm_cov_avail = true;
            else
                IO.spm_cov_avail = false;
            end
            
            %% Check that no temporal filter was used
            if isfield(IO.xX,'K') && isstruct(IO.xX.K)
                mess = GenerateMessageEntry(mess, 'Error in SPM import: No first level analysis with temporal correlations allowed.');
                return
            end

            %% Sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design
            IO.xX = correct_xX(IO.xX);
            
            %% Check whether SPM design matches NM framework
            if numel(SPM.xX.iH) < 2 && strcmp(IO.modeflag,'classification')
                mess = GenerateMessageEntry(mess, 'ERROR in SPM import: You have selected the classification framework but SPM.mat does not contain a group-based design.');
                return
            elseif isempty(SPM.xC) && strcmp(IO.modeflag,'regression')
                mess = GenerateMessageEntry(mess, 'ERROR in SPM import: You have selected the regression framework but SPM.mat does not contain any covariates.');
                return
            end

            %% Check whether SPM contains repeated design
            IO.repeated_anova = ~isempty(IO.xX.iB);
            if IO.repeated_anova
                mess = GenerateMessageEntry(mess, 'ERROR in SPM import: Repeated designs not allowed.');
                return
            end
    
            %% Retrieve mask info
            if isempty(IO.brainmask) || ~exist(IO.brainmask,'file')
                if isfield(SPM,'VM') && isstruct(SPM.VM)
                    IO.Vm = SPM.VM;
                    IO.brainmask = SPM.VM.fname;
                elseif isfield(SPM.xM,'VM') && isstruct(SPM.xM) && ~isempty(SPM.xM.VM)
                    IO.Vm = SPM.xM.VM;
                    IO.brainmask = SPM.xM.VM.fname;
                else
                    IO.Vm = [];
                    IO.brainmask = [];
                end
                if ~isempty(IO.brainmask)
                    IO.Thresh = nk_DataIO3_SpaceDefImage_config(IO.Vm, Thresh.Lm);
                end
            end
    
            %% Generate NM design 
             % Check which columns in the design matrix are dummy variables
            for i=1:size(IO.xX.X,2), uC(i) = numel(unique(IO.xX.X(:,i))); end
            
            switch IO.modeflag
                case 'classification'
                    uC_dummy = uC==2; nuC_dummy = sum(uC_dummy);
                    IO.n_samples = nuC_dummy;
                    if nuC_dummy < 2
                        mess = GenerateMessageEntry(mess, 'ERROR in SPM import: A minimum of 2 dummy variables is required for the classification framework.');
                        return
                    end
                    mess = GenerateMessageEntry(mess, sprintf('* Found %g column(s) with dummy variables in the SPM design.',IO.n_samples), 0, 'blue'); 
                    
                case 'regression'
                    cont = uC>2; n_cont = sum(cont);
                    IO.n_samples=1;
                    if ~n_cont
                        mess = GenerateMessageEntry(mess, 'ERROR in SPM import: No continuous variable found in the SPM design.');
                        return
                    end
                    %IO.n_samples = 1;
                    mess = GenerateMessageEntry(mess, sprintf('Found %g column(s) with continuous variables in the SPM design.',n_cont), 0, 'blue'); 
                    
            end
            
            %% Retrieve global normalization options from SPM
            IO.globscale = SPM.xGX.iGXcalc;
            if ~SPM.xGX.GM
                IO.globnorm = 1;
            else
                IO.globnorm = SPM.xGX.GM;
            end
            switch IO.globscale 
                case 1
                    IO.g=[];
                case {2,3}
                    if isfield(IO,'design')
                        ind = any(IO.design,2);
                    else
                        ind = true(numel(SPM.xGX.gSF),1);   
                    end
                    IO.g = SPM.xGX.gSF(ind);
            end
            [ IO, mess ] = RetrieveImageInfo(IO, IO.datasource, mess);
            IO.n_subjects_all= size(char(IO.PP),1);
            [IO.ID, IO.files, mess] = nk_DefineCaseNames2(IO.PP, sum(IO.n_subjects_all),[], mess);
            
            mess = GenerateMessageEntry(mess, '* SPM design is OK', 0, 'blue');
            
        case 'check_label'
            
            for i=1:size(IO.xX.X,2), uC(i) = numel(unique(IO.xX.X(:,i))); end
            
            switch IO.modeflag
                case 'classification'
                    defsel = IO.xX.iH;
                case 'regression'
                    cont = uC>2;
                    defsel= find(cont);
            end
            
            if ~isfield(IO,'sel_dummy')
                IO.sel_dummy =  nk_input('Specify which columns of the SPM design you want to define as NM labels',0,'e',defsel);  
            end
            
            if isfield(IO,'sel_dummy')
                % Extract columns, files and volumes structures from design matrix
                IO.n_samples = numel(IO.sel_dummy);
                IO.design = SPM.xX.X(:,IO.sel_dummy);
                IO.desc_groups = SPM.xX.name(IO.sel_dummy);
                switch IO.modeflag
                    case 'classification'
                        for i=1:IO.n_samples
                            IO.P{i} = IO.PP(IO.design(:,i)==1,:);
                            IO.V{i} = SPM.xY.VY(IO.design(:,i)==1);
                        end
                    case 'regression'
                        IO.P{1} = char(IO.PP);
                        IO.V{1} = SPM.xY.VY;
                        IO.L = IO.design;
                end
            else
                 IO.n_samples=1;
                 IO.P{1} = char(IO.PP);
                 IO.V{1} = SPM.xY.VY;
                 if isfield(IO,'label')
                    IO.L = IO.label;
                 else
                    IO.L = ones(size(IO.PP,1),1);
                 end
            end
            IO.PP=[]; for i = 1:IO.n_samples, IO.PP = [IO.PP; IO.P{i}]; end
            [ IO, mess ] = RetrieveImageInfo(IO, IO.datasource, mess);
            IO.n_subjects_all= size(char(IO.PP),1);
            [IO.ID, IO.files, mess] = nk_DefineCaseNames2(IO.PP, sum(IO.n_subjects_all),[], mess);
            IO.n_subjects = inf(1,IO.n_samples);
    
        case 'check_cov'
            
                %% Check for covariates 
                if ~isempty(IO.xC) 
                    IO.spm_cov_avail = 1;
                    if ~isfield(IO,'sel_dummy_covs') || isempty(IO.sel_dummy_covs) && isfield(IO,'covmanage_edit') && ~ischar(IO.covmanage_edit) && IO.covmanage_edit < 3
                        nC = numel(SPM.xC);
                        fprintf('\nAvailable covariates in %s:', IO.spm_file);
                        for i=1:nC, fprintf('\n%3g: %s',i, SPM.xC(i).rcname); end
                        IO.sel_dummy_covs =  nk_input(sprintf('Specify which in the SPM design to be defined as covariates in NM(type ''0'' for no covariate extraction, max col number: %g)', nC),0,'e');   
                        if ~(numel(IO.sel_dummy_covs)==1 && ~IO.sel_dummy_covs) && max(IO.sel_dummy_covs)<=nC && numel(IO.sel_dummy_covs) <= nC
                            IO.covars = zeros(size(SPM.xX.X,1),numel(IO.sel_dummy_covs));
                            IO.covnames = cell(1,size(IO.covars,2)); 
                            for i=1:numel(IO.sel_dummy_covs)
                                IO.covars(:,i) = SPM.xC(IO.sel_dummy_covs(i)).c;
                                IO.covnames{i} = SPM.xC(IO.sel_dummy_covs(i)).rcname;
                            end
                        end
                    end
                end   
    end

catch ERR
    mess = GenerateMessageEntry(mess, sprintf('ERROR in SPM import: %s (%s)\nFunction: %s (Line: %g)', ...
        ERR.message, ERR.identifier, ERR.stack(1).name, ERR.stack(1).line));
    if isfield(IO,'sel_dummy'), IO = rmfield(IO,'sel_dummy'); end
    return
end

function xX = correct_xX(xX)

% vector of covariates and nuisance variables
iCG = [xX.iC xX.iG];
iHB = [xX.iH xX.iB];

% set columns with covariates and nuisance variables to zero
X = xX.X;
X(:,iCG) = 0;

ncol = size(X,2);

% calculate sum of columns
% The idea behind this is that for each factor the sum of all of its columns should be "1".
Xsum = zeros(size(X));
for i=1:ncol
  % only sum up columns without covariates and nuisance variables
  if isempty(find(iCG==i))
    Xsum(:,i) = sum(X(:,1:i),2);
  end
end

% find columns where all entries are constant except zeros entries
% that indicate columns with covariates and nuisance variables
ind = find(any(diff(Xsum))==0 & sum(Xsum)>0);

% no more than 2 factors expected
if length(ind) > 2
  error('Weird design was found that cannot be analyzed correctly.');
end

% correction is only necessary if 2 factors (iH/iB) were found
if length(ind) > 1
  iF = cell(length(ind),1);

  j = 1;
  % skip columns with covariates and nuisance variables
  while find(iCG==j),  j = j + 1; end

  for i=j:length(ind)
    iF{i} = [j:ind(i)];
  
    j = ind(i)+1;
    % skip columns with covariates and nuisance variables
    while find(iCG==j), j = j + 1; end
  end
  
  % not sure whether this will always work but usually iB (subject effects) should be longer than iH (time effects)
  if length(iF{1}) > length(iF{2})
    xX.iB = iF{1};
    xX.iH = iF{2};
  else
    xX.iB = iF{2};
    xX.iH = iF{1};
  end
end