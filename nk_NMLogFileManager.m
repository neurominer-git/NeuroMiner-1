function [status, analysis] = nk_NMLogFileManager(act, dat, analysis, caller, mess)

myname = 'NM Logfile Manager';
if ~exist('caller','var'), caller = []; end
if ~exist('mess','var'), mess = []; end
if strcmp(act,'add_entry') && ( isempty(caller) || isempty(mess) )
     status.error.message = 'Logfile manager cannot add an empty message to the logfile';
     status.error.identifier = 'NM:logfilemanager:emptymessage';
     status.error.function = 'nk_NMLogFileManager';
     status.error.line = 6;
     return
end
status = [];

try
    switch act
        case 'check_exist'
            if isfield(analysis,'logfile') && exist(analysis.logfile,'file'), status = true; end
        case 'init'
            wrtfl = 'w';
            analysis.logfile = fullfile(analysis.rootdir,sprintf('NM_Analysis_%s.log',analysis.id));
            fl = nk_NMLogFileManager('check_exist', dat, analysis);
            analysis.paramdir = fullfile(analysis.rootdir,'params');
            analysis.paramfile = fullfile(analysis.paramdir,sprintf('NM_Analysis_%s_params_%s.mat',analysis.id,timestampstr));
            if ~exist(analysis.paramdir,'dir'), mkdir(analysis.paramdir); end
            params = analysis.params; meta = analysis.meta; id = analysis.id;
            save(analysis.paramfile,'params', 'meta', 'id')
            existfl=false;
            if islogical(fl) && fl
                %prompt = {'NM found a previous analysis log file!',sprintf('%s',analysis.logfile)};
                %exist_reply = questdlg(prompt,myname,'Overwrite','Append','Abort','Abort');
%                 switch exist_reply
%                     case 'Append'
%                         wrtfl = 'a';
%                     case 'Abort'
%                         return
%                 end
                wrtfl = 'a';
                existfl = true;
            elseif isstruct(fl)
                errormsg = {status.error.message, status.error.identifier};
                errordlg(errormsg,'NM Logfile Manager');
            end
            fid = fopen(analysis.logfile,wrtfl);
            if existfl,
                fprintf(fid,'\nOverwrite previous analysis settings!\n');
            end
            fprintf(fid,'NM ANALYSIS LOG FILE');
            fprintf(fid,'\n====================');
            fprintf(fid,'\nInit date:\t\t\t%s',             analysis.meta.TIME);
            fprintf(fid,'\nBy user:\t\t\t%s',               char(analysis.meta.USER));
            fprintf(fid,'\nOS Type:\t\t\t%s',               char(analysis.meta.OS.name));
            fprintf(fid,'\nOS Version:\t\t\t%s',            char(analysis.meta.OS.version));
            fprintf(fid,'\nOS Architecture:\t%s',           char(analysis.meta.OS.arch));
            fprintf(fid,'\nMATLAB version:\t\t%s',          version);
            fprintf(fid,'\nNM version:\t\t\t%s',            analysis.meta.NM.ver);
            fprintf(fid,'\nNM workspace ID:\t%s',           analysis.params.id);
            fprintf(fid,'\nAnalysis ID:\t\t%s',             analysis.id);
            fprintf(fid,'\nAnalysis root dir:\t%s',         analysis.rootdir);
            fprintf(fid,'\nAnalysis description:');
            for i=1:numel(analysis.desc)
                fprintf(fid,'\n\t\t\t\t\t%s',analysis.desc{i});
            end
            fprintf(fid,'\n\nAnalysis overview:');
            fprintf(fid,'\n==================\n');
            str = nk_GetAnalysisInfo(dat, analysis);
            for i=1:numel(str), fprintf(fid,'%s',str{i}); end
            fclose(fid);
        case 'delete'
            fl = nk_NMLogFileManager('check_exist', dat, analysis);
            if islogical(fl) && fl
                prompt = {'Are you sure you want to delete the analysis logfile:',sprintf('%s',analysis.logfile)};
                delete_reply = questdlg(prompt,myname,'Yes','No');
                if strcmp(delete_reply,'Yes'), delete(analysis.logfile); end
            end
        case 'add_entry'
            fid = fopen(analysis.logfile,'a');
            fprintf(fid,'\n\nLogfile entry time: %s',datestr(now));
            fprintf(fid,'\n=========================================');
            if isobject(mess)
                fprintf(fid,'\n* Logfile entry type: ERROR');
                fprintf(fid,'\n* Error message:\t\t%s',mess.message);
                fprintf(fid,'\n* Error identify:\t\t%s',mess.identifier);
                fprintf(fid,'\n* Error stack:');
                for i=1:numel(mess.stack)
                    fprintf(fid,'\n\t in %s, name: %s, line: %g.', ...
                        mess.stack(i).file,mess.stack(i).name,mess.stack(i).line);
                end
            elseif isstruct(mess) && isfield(mess,'text')
                fprintf(fid,'\n%s',mess.text);
            end     
    end
catch ERR
    status = [];
    status.error.message = ERR.message;
    status.error.identifier = ERR.identifier;
    status.error.function = ERR.stack(1).name;
    status.error.line = ERR.stack(1).line;
end