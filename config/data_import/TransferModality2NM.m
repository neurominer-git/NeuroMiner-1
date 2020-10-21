function [D, I] = TransferModality2NM(D, I, varind)

% Write IO data to D
if (~isfield(I,'oocvflag') || ~I.oocvflag) && isfield(I,'modeflag'), D.modeflag = I.modeflag; end

if isfield(I,'Y')
    
    % Check whether multiple modalities have to be imported
    if iscell(I.Y)
        nM = numel(I.Y)-1; Y = I.Y;
        m = size(I.Y{1},1);
    else
        nM = 0; Y{1} = I.Y; 
        m = size(I.Y,1);
    end
    cnt = 1;
    
    % Select cases based on index vector 'selCases'
    if ~isfield(I,'selCases') || isempty(I.selCases) || numel(I.selCases)~=numel(I.ID)
        selCases = true(m,1);
    else
        selCases = I.selCases;
    end
    I.ID = I.ID(selCases,:);
    
    % Perform integrity check of cases
    if isfield(D,'cases')
        flg = nk_DefineCaseListIntegrity(D.cases, I.ID) ;
        if ~flg, return; end
    end
    
    % Loop through modalities to be imported
    for i=varind:nM+varind
        
        % Map current modality data to iY
        if iscell(I.Y),
            n = size(I.Y{cnt},2);
            iY = I.Y{cnt};
            if isfield(I,'Yw'), iYw = I.Yw{cnt}; else, iYw=[];end
        else
            n = size(I.Y,2);
            iY = I.Y;
            if isfield(I,'Yw'), iYw = I.Yw; else, iYw=[]; end
        end
        
        % Select features accordin to index vector 'selFeats'
        if isfield(I,'selFeats') && ~isempty(I.selFeats)
            selFeats = I.selFeats;
        else
            selFeats = true(1,n);
        end
        
        % Retrieve feature names, respectively 
        if isfield(I,'featnames') && ~isempty(I.featnames)
            featnames = I.featnames(selFeats);
        else
            featnames = [];
        end
        
        % Noe extract the data for current modality accordingly
        D.Y{i} = iY(selCases, selFeats);
        
        % Save further modality-specific information
        switch I.datasource
            case {'surf','spm','nifti'}
                if strcmp(I.datasource,'surf')
                    D.datadescriptor{i}.type = 2;
                else
                    D.datadescriptor{i}.type = 1;
                end
                D.datadescriptor{i}.source = 'image';
                D.datadescriptor{i}.globscale = I.globscale;
                D.datadescriptor{i}.globnorm = I.globnorm;
                if isfield(I,'g'), D.datadescriptor{i}.g = I.g; end;
                D.datadescriptor{i}.threshop = I.Thresh.threshop;
                D.datadescriptor{i}.threshval = I.Thresh.Vml(cnt);
                if ~isempty(iYw), D.datadescriptor{i}.Yw = iYw(:,selFeats); end
                if isfield(I.Thresh,'Lm') && numel(I.Thresh.Lm)>1
                    D.datadescriptor{i}.desc = sprintf('%s [ %s ]', I.desc_data, I.Thresh.Lm{cnt}) ;
                    D.datadescriptor{i}.desc_label = I.Thresh.Lm{cnt};
                else
                    D.datadescriptor{i}.desc = I.desc_data;
                    D.datadescriptor{i}.desc_label = [];
                end
                D.featnames{i} = [];
                D.datadescriptor{i}.Vmvol = I.Vmvol;        
                D.datadescriptor{i}.Vm = I.Vm;
            case {'matrix','vector'}
                D.datadescriptor{i}.type = 0;
                D.datadescriptor{i}.source = 'matrix';
                D.datadescriptor{i}.desc = I.desc_data;
                D.featnames{i} = featnames;
        end
        if isfield(I,'brainmask'), D.brainmask{i} = I.brainmask; end;
        D.badcoords{i} = false(1,size(D.Y{i},2)); 
        D.datadescriptor{i}.input_settings = I;
        D.datadescriptor{i}.input_settings = rmfield(D.datadescriptor{i}.input_settings,'Y');
        if isfield(I,'PP'), D.files{i} = I.PP; else, D.files{i} = []; end
        cnt = cnt+1;
    end
    
    D.cases = I.ID;
    if I.labels_known
        
        if ~isfield(D,'groupnames'), D.groupnames = I.desc_groups; end
        if isfield(I,'desc_labels'), D.labelnames = I.desc_labels; end
        if ~isfield(D,'n_subjects')
            if iscell(I.n_subjects)
                if any(cellfun(@any,cellfun(@isfinite,I.n_subjects,'UniformOutput',false))), D.n_subjects = I.n_subjects; end
            else
                if ~any(~isfinite(I.n_subjects)), D.n_subjects = I.n_subjects; end
            end     
        end
        if ~isfield(D,'n_subjects_all') || ~isfinite(D.n_subjects_all), D.n_subjects_all = I.n_subjects_all; end
        if ~isfield(D,'label'), 
            if isfield(I,'label')
                if numel(I.label)==numel(selCases)
                   D.label = I.label(selCases); 
                else
                   D.label = I.label;
                end
            else
                D.label = I.L;
            end
        elseif isfield(I,'labelmanage_edit') && I.labelmanage_edit == 1
            D.label = [D.label I.label]; 
        end
    end
    
    if isfield(I,'covars')
        if isfield(D,'covars')
            switch I.covmanage_edit
               case 1
                   D.covars = [D.covars I.covars];
                   D.covnames = [D.covnames I.covnames];
               case 2
                   D.covars = I.covars;
                   D.covnames = I.covnames;
            end
        else
            try
                D.covars = I.covars;
                D.covnames = I.covnames;
            catch
                warning('Covnames missing')
            end
       end
    end
    
    if isfield(I,'survanal_time') && isnumeric(I.survanal_time)
        D.time = I.survanal_time;
    end
    
    if ~isfield(D,'defs') || ~isfield(D.defs,'os_sys')
        D.defs.os_sys = computer;
        if isunix
            [~,D.defs.os_ver] = unix('sw_vers');
        elseif ispc
            D.defs.os_ver = system_dependent('getwinsys');
        end
        D.defs.matlab_ver = version;
    end
    if nargout==2
        I = rmfield(I,'Y');
        if isfield(I,'s_ID'); I = rmfield(I,'s_ID'); end
        if isfield(I,'s_label'); I = rmfield(I,'s_label'); end
        if isfield(I,'label'); I = rmfield(I,'label'); end
        if isfield(I,'t_Y'); I = rmfield(I,'t_Y'); end
        if isfield(I,'s_featnames'); I = rmfield(I,'s_featnames'); end
        if isfield(I,'Yw'), I=rmfield(I,'Yw'); end
    end
end