% UIWAIT makes nk_PrintResults2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function load_popupmenu1(handles)

switch handles.modeflag
    case 'classification'
        for i = 1:handles.nclass
            popuplist{i} = sprintf('Classifier %g: %s', i, handles.BinClass{i}.description);
        end
        if isfield(handles,'MultiClass'), popuplist{end+1} = sprintf('Multi-group classifier'); end
    case 'regression'
         popuplist{1} = sprintf('Regression model');
end

handles.popupmenu1.String = popuplist;
