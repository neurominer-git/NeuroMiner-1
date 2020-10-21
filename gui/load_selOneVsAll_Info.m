function load_selOneVsAll_Info(handles)

popuplist{1} = sprintf('Multi-group display');
for i = 1:numel(handles.NM.groupnames)
    popuplist{end+1} = sprintf('%s vs. REST', handles.NM.groupnames{i});
end
set(handles.selOneVsAll_Info, 'String', popuplist);