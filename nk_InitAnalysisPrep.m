function [act, NM, A] = nk_InitAnalysisPrep(NM, A, parentstr)

disallow = false;
analdim = 1;
amode = 1;
ovrwrt = 3;
na_str = '?'; mn_str = []; mn_act = [];
Amodes = {'generate new','manage existing'};
Aovrwrt = {'delete', 'delete & wipe', 'complete reset', 'parameter reset', 'update paths and descriptors', 'overwrite NM parameter template'};

if ~exist('A','var') || isempty(A)
    A.mode = amode;
    A.analdim = analdim;
    A.ovrwrt = ovrwrt;
    A.id = na_str;
    A.desc = na_str;
    A.parentdir = na_str;
    A.rootdir = na_str;
end

if isfield(NM,'analysis')
    
    mode_str = Amodes{A.mode}; 

    if A.mode == 2
        if A.analdim>numel(NM.analysis), A.analdim=numel(NM.analysis); end
        mn_sel_mode = sprintf('Generate new or manage existing analyses [ %s ]|', mode_str );
        mn_str = mn_sel_mode; mn_act = {'manage_analyses'}; 
        if isfield(NM.analysis{A.analdim},'id')
            if strcmp(A.id,na_str), A.id = NM.analysis{A.analdim}.id; end
            if strcmp(A.desc,na_str), A.desc = NM.analysis{A.analdim}.desc; end
            if strcmp(A.parentdir,na_str), A.parentdir = NM.analysis{A.analdim}.parentdir; end
            if strcmp(A.rootdir,na_str), A.rootdir = NM.analysis{A.analdim}.rootdir; end
        end
        anal_str = sprintf('Analysis %g',A.analdim);
        mn_sel_anal = sprintf('Select existing analysis [ %s ]|', anal_str);
        mn_str = [mn_str mn_sel_anal]; mn_act = [mn_act 'sel_analyses']; 
        
        todo_str = Aovrwrt{A.ovrwrt};
        mn_todo_anal = sprintf('Specify what to do with selected analysis [ %s ]|', todo_str);
        mn_str = [mn_str mn_todo_anal]; mn_act = [mn_act 'todo_analyses']; 
    else
        A.analdim = numel(NM.analysis)+1;
        mn_sel_mode = sprintf('Generate new or manage existing analyses [ %s: analysis %g ]|', mode_str, A.analdim );
        mn_str = mn_sel_mode; mn_act = {'manage_analyses'}; 
    end
end

if A.ovrwrt > 1 && A.ovrwrt<6
    if ~strcmp(A.id,na_str), def_str = A.id; else, def_str = na_str; end
    mn_def_analysis_id = sprintf('Define analysis identifier [ %s ]|', def_str);
    if isempty(mn_act)
        mn_str = mn_def_analysis_id; mn_act = {'def_analysis_id'};
    else
        mn_str = [mn_str mn_def_analysis_id]; mn_act = [mn_act 'def_analysis_id']; 
    end

    if ~strcmp(A.desc,na_str), desc_str = 'provided'; else, desc_str = na_str; end
    mn_def_analysis_desc = sprintf('Provide analysis description [ %s ]|', desc_str);
    mn_str = [mn_str mn_def_analysis_desc]; mn_act = [mn_act 'def_analysis_desc']; 
    
    mn_def_parentdir = sprintf('Specify parent directory of the analysis [ %s ]|', spm_file(A.parentdir,'short30'));
    mn_str = [mn_str mn_def_parentdir]; mn_act = [mn_act 'def_parentdir']; 
    
end

if any(strcmp(A.parentdir,na_str)) || any(strcmp(A.desc,na_str)) || any(strcmp(A.id,na_str)), 
    disallow=true; 
end

if ~disallow
    mn_run = sprintf('PROCEED >>>'); mn_str = [ mn_str mn_run ] ; mn_act = [mn_act 'run'];
end
nk_PrintLogo
fprintf('\n'); mestr = 'Select action of the analysis manager';  navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); cprintf('*blue','You are here: %s >>> ',parentstr); 
act = char(nk_input(mestr,0,'mq', mn_str, mn_act));

switch act
    case 'BACK'
        return
    
    case 'manage_analyses'
        modus = nk_input('Define analysis manager operation mode',0,'mq','generate new|manage existing',[1,2],A.mode);
        if modus>0, A.mode = modus; end
    
    case 'sel_analyses'
        analdim = A.analdim;
        t_act = 1; brief = 1; while t_act>0, [t_act, analdim, NM, ~, brief ]= nk_SelectAnalysis(NM, 0, 'MAIN >>> Select Analysis', analdim, 0,0,[],brief); end
        if ~isempty(analdim)
            A.analdim = analdim;
            if numel(analdim)==1
                A.id = NM.analysis{A.analdim}.id ;
                A.desc = NM.analysis{A.analdim}.desc;             
                A.parentdir = NM.analysis{A.analdim}.parentdir;
            end
        end
            
    case 'todo_analyses'
        ovrwrt = nk_input('Define what to do with the selected analysis',0,'mq', ...
                        ['delete|' ...
                         'delete and wipe from computer|' ...
                         'completely reset|' ...
                         'reset parameters (risk of inconsistency)|' ...
                         'update paths and descriptors|' ...
                         'overwrite current NM parameter template with analysis workspace'],1:6,A.ovrwrt);
        if ovrwrt>0, A.ovrwrt = ovrwrt; end
        
    case 'def_parentdir'
        
         parentdir = uigetdir(pwd,'Select parent directory of the new analysis');
         if ischar(parentdir) && exist(parentdir,'dir')==7, A.parentdir = parentdir; end
             
    case 'def_analysis_id'
        if isfield(A,'id') && ~isempty(A.id), id_def = A.id;  end
        id = nk_input('Provide an alphanumeric identifier for your analysis',0,'s',id_def);
        if isempty(id) || strcmp(id,''), A.id=na_str; return; end
        A.id = id;
        ind_alphanum = isstrprop(A.id,'alphanum');
        ind_minus = false(size(ind_alphanum)); ind_dot=ind_minus;ind_under=ind_minus;
        f_minus = strfind(A.id,'-'); if ~isempty(f_minus);ind_minus(f_minus)=true; end
        f_under = strfind(A.id,'_'); if ~isempty(f_under);ind_under(f_under)=true; end
        f_dot = strfind(A.id,'.');   if ~isempty(f_dot);ind_dot(f_dot)=true; end
        ind = ~(ind_alphanum | ind_minus | ind_dot | ind_under);
        A.id(ind) = '';
        
    case 'def_analysis_desc'
        if isfield(A,'desc') && ~isempty(A.desc), desc_def = A.desc;  end
        A.desc = nk_input('Describe the analysis you want to initialize',0,'s+',desc_def);
        
    case 'run'
        
        switch A.ovrwrt
    
            case {1,2}
                if A.ovrwrt == 2
                    askfl = questdlg('Are you sure you want to wipe this analysis from your computer?',mestr,'Yes','No','No');
                    if strcmp(askfl,'Yes'), 
                        for i=1:numel(A.analdim)
                            rmdir( NM.analysis{A.analdim(i)}.rootdir,'s' ); 
                        end
                    end
                end
                NM.analysis(A.analdim) = [];
                if isempty(NM.analysis), NM = rmfield(NM,'analysis'); end
            case {3,4}
                if A.ovrwrt == 3
                    NM.analysis{A.analdim} = []; 
                    NM.analysis{A.analdim}.id                 = A.id;                   
                    NM.analysis{A.analdim}.desc               = A.desc;
                    NM.analysis{A.analdim}.parentdir          = A.parentdir;
                    NM.analysis{A.analdim}.rootdir            = fullfile(NM.analysis{A.analdim}.parentdir, sprintf('NM_ID%s_A%g_%s',NM.id,A.analdim,NM.analysis{A.analdim}.id));
                    if ~exist(NM.analysis{A.analdim}.rootdir,'dir'), mkdir(NM.analysis{A.analdim}.rootdir); end
                end
                NM.analysis{A.analdim}.params.TrainParam      = NM.TrainParam;
                NM.analysis{A.analdim}.params.datadescriptor  = NM.datadescriptor;
                NM.analysis{A.analdim}.params.modeflag        = NM.modeflag;
                NM.analysis{A.analdim}.params.cv              = NM.cv;
                NM.analysis{A.analdim}.params.id              = NM.id;
                NM.analysis{A.analdim}.meta.TIME              = datestr(now);
                NM.analysis{A.analdim}.meta.USER              = java.lang.System.getProperty('user.name');
                NM.analysis{A.analdim}.meta.OS.name           = java.lang.System.getProperty('os.name');
                NM.analysis{A.analdim}.meta.OS.version        = java.lang.System.getProperty('os.version');
                NM.analysis{A.analdim}.meta.OS.arch           = java.lang.System.getProperty('os.arch');
                NM.analysis{A.analdim}.meta.MATLAB.ver        = ver;
                NM.analysis{A.analdim}.meta.NM.ver            = NM.defs.NM_ver;
                NM.analysis{A.analdim}.status = 0;
                [ log_status, NM.analysis{A.analdim} ]        = nk_NMLogFileManager('init', NM, NM.analysis{A.analdim});
            case 5
                old_rootdir = sprintf('NM_ID%s_A%g_%s',NM.analysis{A.analdim}.params.id,A.analdim,NM.analysis{A.analdim}.id); 
                old_rootdirpath = NM.analysis{A.analdim}.rootdir;
                new_rootdirpath = fullfile(A.parentdir,old_rootdir);
                if  exist(new_rootdirpath,'dir') && ~strcmp(old_rootdirpath, new_rootdirpath)
                    NM.analysis{A.analdim}.desc                   = A.desc;
                    NM.analysis{A.analdim}.parentdir              = A.parentdir;
                    NM.analysis{A.analdim}.rootdir                = new_rootdirpath;
                    NM.analysis{A.analdim}.paramdir               = fullfile(new_rootdirpath,'params');
                    [~,n,e]                                       = fileparts(NM.analysis{A.analdim}.paramfile);
                    NM.analysis{A.analdim}.paramfile              = fullfile(NM.analysis{A.analdim}.paramdir,n,e);
                    NM.analysis{A.analdim}.logfile                = fullfile(new_rootdirpath,sprintf('NM_Analysis_%s.log', NM.analysis{A.analdim}.id));
                    mess.text = sprintf('Description:');
                    for i=1:numel(A.desc)
                        mess.text = sprintf('%s\n\t\t\t\t\t%s', mess.text, A.desc{i});
                    end
                    mess.text = sprintf('%s\nParent directory: %s\nRoot directory: %s', mess.text, A.parentdir, new_rootdirpath);
                    log_status = nk_NMLogFileManager('add_entry', NM, NM.analysis{A.analdim}, 'InitAnalysis:UpdatePaths', mess);
                end
            case 6
                NM.TrainParam = NM.analysis{A.analdim}.params.TrainParam;
                NM.cv = NM.analysis{A.analdim}.params.cv;
        end
        act='BACK';
end
       


