function Status = nk_GetAnalysisStatus(NM)

Status.completed_analyses  = []; 
Status.isequal_cv          = []; 
Status.nmodal_analyses     = [];
Status.analexistflag       = false; 
Status.analreadyflag       = false; 
Status.analcompleteflag    = false;
Status.oocvreadyflag       = false;
Status.oocvappflag         = false;

if isfield(NM,'analysis')
    n_anal = numel(NM.analysis);
    Status.completed_analyses = false(1,n_anal);
    Status.isequal_cv = false(1,n_anal);
    Status.nmodal_analyses = zeros(1,n_anal);
    for i = 1:n_anal
        if NM.analysis{i}.status, 
            Status.completed_analyses(i)= true; 
        end
        if isequaln(NM.cv, NM.analysis{i}.params.cv),
            Status.isequal_cv(i) = true; 
        end
        if Status.completed_analyses(i)
            Status.nmodal_analyses(i) = numel(NM.analysis{i}.GDdims);
        end
    end
end

if ~isempty(Status.completed_analyses),    
    Status.analexistflag = true; 
    if any(Status.completed_analyses), Status.analreadyflag = true; end    
    if sum(Status.completed_analyses) == numel(Status.completed_analyses), Status.analcompleteflag = true; end
    if Status.analcompleteflag && isfield(NM,'OOCV'), Status.oocvreadyflag = true; end
    if Status.oocvreadyflag && isfield(NM.defs,'analyses_locked') && NM.defs.analyses_locked, Status.oocvappflag = true; end
end

