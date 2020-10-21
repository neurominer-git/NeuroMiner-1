function load_selCVoocv(handles)

popupstr{1} = 'Training/CV data';
fnd=false;
if isfield(handles,'OOCVinfo') && isfield(handles.OOCVinfo,'AnalVec') && ismember(handles.curranal, handles.OOCVinfo.AnalVec)
    handles.selCVoocv.Enable = 'on';
    fI = find(handles.OOCVinfo.AnalVec==handles.curranal); cnt=2;
    if handles.OOCVinfo.Analyses{fI}.OOCVdone
        for i = 1:numel(handles.OOCVinfo.Analyses{fI}.OOCVvec)
            popupstr{cnt} = handles.OOCVinfo.Analyses{fI}.descriptor{i};
            cnt=cnt+1; fnd=true;
        end
    end
    if ~fnd,  handles.selCVoocv.Enable = 'off'; end
else
    handles.selCVoocv.Enable = 'off';
end
handles.selCVoocv.String = popupstr;