function set_panel_visibility(handles,flag)

% Visualization panel
set(handles.pnVisual,'Visible',flag);
set(handles.pn3DView,'Visible',flag);
% Classification / Regression plot
set(handles.pnBinary,'Visible',flag);
% Model performance
set(handles.pnModelPerf,'Visible',flag);
% Binary classifier selector
set(handles.popupmenu1,'Enable',flag);
% Prediction metric selector
set(handles.selModelMeasures,'Enable',flag);

