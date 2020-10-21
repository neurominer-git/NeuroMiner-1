function load_selSubParams(handles)

curclass = get(handles.popupmenu1,'Value');
if curclass>numel(handles.ModelParamsDesc), curclass=1; end
popupstr = handles.ModelParamsDesc{curclass}';
popupstrgr1{1} = 'All parameters'; cnt=2;
for i=1:numel(popupstr)
    if iscell(handles.ModelParams{curclass})
        nuPi = numel(unique(cell2mat(handles.ModelParams{curclass}(:,i))));
    else
        nuPi = numel(unique(handles.ModelParams{curclass}(:,i)));
    end
    if nuPi>1
        popupstrgr1{cnt}=popupstr{i};
        cnt=cnt+1;
    end
end
handles.selSubParam.String = popupstrgr1' ;