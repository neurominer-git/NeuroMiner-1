function handles = load_modal(handles, GDdims)

for i=1:numel(GDdims)
    popuplist{i} = sprintf('Modality #%g: %s',i,GDdims{i}.datadescriptor.desc);
end
set(handles.selModal, 'String', popuplist);