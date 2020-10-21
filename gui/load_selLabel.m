function load_selLabel(handles)
    
labelstr='';
for i=1:size(handles.NM.label,2)
    if isfield(handles.NM,'labelnames')
        labelstr = handles.NM.labelnames{i};
    end
      popuplist{i} = sprintf('Label #%g%s',i,labelstr);
end
set(handles.selLabel, 'String', popuplist);