function STATUS = nk_CheckFieldStatus(parent, children, grandchildren, grandchildrencrit, STATUS)

strdef = '...'; strundef = '???'; strcrit = '!!!';

n = numel(grandchildren);
missdef = cellstr(repmat(strundef,n,1));
okdef   = cellstr(repmat(strdef,n,1));
if exist('grandchildren','var') && ~isempty(grandchildren)
    if exist('grandchildrencrit','var') && ~isempty(grandchildrencrit)
        missdef{grandchildrencrit} = strcrit;
    end
else
    grandchildren = [];
end

if isempty(parent)
    childrenfields = [];
    nC = 0;
else
    childrenfields = fieldnames(parent);
    nC = numel(childrenfields);
end

newfl = 1;
for i=1:numel(children)
    STATUS.(children{i}) = strundef;
    for j=1:nC
        if strcmp(childrenfields{j},children{i}) 
            newfl = 0; STATUS.(children{i}) = strdef; break
        end
    end
end

if ~isempty(grandchildren) && nC >0
    grandchildrenfields = fieldnames(parent.(children{1}));
    for i=1:numel(grandchildren)
        if newfl
            STATUS.(grandchildren{i}) = missdef{i};
        else
            fl = 1;
            for j=1:numel(grandchildrenfields)
                if strcmp(grandchildrenfields{j},grandchildren{i}); fl = 0; break; end
            end
            if fl
                STATUS.(grandchildren{i}) = missdef{i};
            else
                STATUS.(grandchildren{i}) = okdef{i};
            end
        end
    end
end
%     if newfl || ~isfield(NM,'cv') || isempty(NM.cv), STATUS.CV = 'CV'; else STATUS.CV= '...'; end
%     if newfl || ~isfield(NM,'SVM') || isempty(NM.SVM), STATUS.SVM = 'SVM'; else STATUS.SVM = '...'; end
%     if newfl || ~isfield(NM.TrainParam,'GRD') || isempty(NM.TrainParam.GRD), STATUS.GRD = 'UNDEFINED'; else STATUS.GRD = '...'; end
%     if newfl || ~isfield(NM.TrainParam,'PREPROC') || isempty(NM.TrainParam.PREPROC), 
%         STATUS.PREPROC = 'UNDEFINED'; STATUS.FEATGEN = 'UNDEFINED'; 
%     else
%         STATUS.PREPROC = '...'; STATUS.FEATGEN = '...'; 
%     end
%     if newfl || ~isfield(NM.TrainParam,'RFE') || isempty(NM.TrainParam.RFE), STATUS.RFE = 'UNDEFINED'; else STATUS.RFE = '...'; end
%     if newfl || ~isfield(NM.TrainParam,'MULTI') || isempty(NM.TrainParam.MULTI), STATUS.MULTI = 'UNDEFINED'; else STATUS.MULTI = '...'; end
%     if newfl || ~isfield(NM.TrainParam,'VIS') || isempty(NM.TrainParam.VIS), STATUS.VIS= 'UNDEFINED'; else STATUS.VIS = '...'; end
%     if newfl || ~isfield(NM.TrainParam,'SAV') || isempty(NM.TrainParam.SAV), STATUS.SAV = 'UNDEFINED'; else STATUS.SAV = '...'; end
%     if newfl || ~isfield(NM.TrainParam,'OOCV') || isempty(NM.TrainParam.OOCV), STATUS.OOCV = 'UNDEFINED'; else STATUS.OOCV = '...'; end
