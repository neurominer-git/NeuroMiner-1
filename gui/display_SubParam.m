function [ handles, contfl] = display_SubParam(handles, caller)

curclass = handles.popupmenu1.Value;
if strcmp(handles.popupmenu1.String{curclass},'Multi-group classifier'), curclass=1; end
param = handles.selSubParam.String{handles.selSubParam.Value};
Pind = strcmp(handles.ModelParamsDesc{curclass},param);
contfl = true;
if ~any(Pind)
    if exist('caller','var') && ~isempty(caller)
        return
    else
        handles = display_modelperf(handles);
    end
else
    if size(handles.currmeas,2)>1,
        [mPerf, sdPerf,~,bars]=extract_subparam_performance(handles.ModelParams{curclass}, ...
                            handles.currmeas, Pind, [], 1, handles.axes17);
    else
        [mPerf, sdPerf,~,bars]=extract_subparam_performance(handles.ModelParams, ...
                            handles.currmeas, Pind, curclass, 1, handles.axes17);
    end
    handles.axes17.XLabel.String = sprintf('Parameter subspace selection: %s', param);
    if size(handles.currmeas)>1
        bars(1).FaceColor='b';
        bars(2).FaceColor='r';
        handles.legend_modelperf.String =  {'CV1 performance','CV2 performance'};
        handles.legend_modelperf.Visible='on';
    else
        bars(1).FaceColor=rgb('green');
        handles.legend_modelperf.Visible = 'off';
    end
    contfl=false;
    assignin('base', 'nm_viewer_data_mean', mPerf);
    assignin('base', 'nm_viewer_data_sd', sdPerf);
end