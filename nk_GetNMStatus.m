function [status, paramstr] = nk_GetNMStatus(NM)

status = struct(  'import_finished', false, ...
                    'setup_ok', false, ...
                    'analyses_exist', false, ... 
                    'analyses_ready', false, ... 
                    'analyses_completed', false, ...
                    'analyses_locked', false, ...
                    'oocv_data_ready', false, ...
                    'oocv_anal_ready', false);

paramstr=[];
if isfield(NM.defs,'import_finished') && NM.defs.import_finished
   status.import_finished = true;
end
if isfield(NM,'Y')
    [pstatus, paramstr] = nk_SetupGlobVars2(NM,'check');
    if ~pstatus
       status.setup_ok = true;
    end
    AS = nk_GetAnalysisStatus(NM);
    if AS.analexistflag, status.analyses_exist = true; end
    if AS.analreadyflag, status.analyses_ready = true; end
    if AS.analcompleteflag, 
        status.analyses_completed = true; 
    end
    status.completed_analyses = AS.completed_analyses; 
    status.isequalcv_analyses = AS.isequal_cv;
    status.nmodal_analyses = AS.nmodal_analyses;
    if isfield(NM.defs,'analyses_locked') && NM.defs.analyses_locked
        status.analyses_locked = true;
    end
    if AS.oocvreadyflag, 
        status.oocv_data_ready = true;
        if isfield(NM.TrainParam,'OOCV'), status.oocv_anal_ready = true; end
    end
else
    status.import_finished = false;
end
    