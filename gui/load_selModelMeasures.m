function load_selModelMeasures(handles)

l=1;
switch handles.modeflag
    case 'classification'
        popupstr{l} = 'Classification plot'; l=l+1;
    case 'regression'
        popupstr{l} = 'Regression plot'; l=l+1;
end

predind = get(handles.popupmenu1,'Value');
predstr = get(handles.popupmenu1,'String');

if strcmp(handles.METAstr,'none') && isfield(handles,'grid') && sum(handles.MLparams.NumParamCombs) > handles.nclass
    GridVarNames = fieldnames(handles.grid);
    if strcmpi(predstr{predind},'Multi-group classifier')
        for i=1:numel(GridVarNames)
             switch GridVarNames{i}
                    case 'MultiCVPerf'
                        popupstr{l} = 'Multi-class cross-validation performance'; l=l+1;
                    case 'MultiERR_CVTSPerf'
                        popupstr{l} = 'Multi-class generalization error'; l=l+1;
                    case 'MultiComplexity'
                        popupstr{l} = 'Multi-class complexity'; l=l+1;
                    case 'MultiCVDiversity'
                        popupstr{l} = 'Multi-class ensemble diversity'; l=l+1;
                    case 'MultiSelNodeFreq'
                        popupstr{l} = 'Multi-class model selection frequency'; l=l+1;
             end
        end
    else
        for i=1:numel(GridVarNames)
            switch GridVarNames{i}
                case 'mean_CVPerf'
                    popupstr{l} = 'Cross-validation performances'; l=l+1;
                case 'mean_Err_CVTSPerf'
                    popupstr{l} = 'Generalization error'; l=l+1;
                case 'mean_Complexity'
                    popupstr{l} = 'Model complexity'; l=l+1;
                case 'mean_CVDiversity'
                    popupstr{l} = 'Ensemble diversity'; l=l+1;
                case 'SelNodeFreq'
                    popupstr{l} = 'Model selection frequency'; l=l+1;
                case 'mean_SeqGain'
                    popupstr{l} = 'Overall sequence gain'; l=l+1;
                case 'mean_SeqExamFreq'
                    popupstr{l} = 'Examination frequencies'; l=l+1;
                case 'mean_SeqPercUpper'
                    popupstr{l} = 'Case propagation thresholds'; l=l+1;
            end
        end
    end
    handles.selModelMeasures.Enable = 'on';
elseif strcmp(handles.METAstr,'bagged')
    if strcmpi(predstr{predind},'Multi-group classifier')
        popupstr = {popupstr{1},'Multi-class cross-validation performance summary','Multi-class generalization error summary'};
    else
        popupstr = {popupstr{1},'Cross-validation performance summary','Generalization error summary', 'Model complexity summary'};
    end
    handles.selModelMeasures.Enable = 'on';
else
    handles.selModelMeasures.Enable = 'off';
end

if isfield(handles,'visdata') && ~isempty(handles.visdata) && strcmp(handles.METAstr,'none')
    popupstr{l} = 'Visualization results'; 
    handles.selModelMeasures.Enable = 'on';
end

handles.selModelMeasures.String = popupstr;
