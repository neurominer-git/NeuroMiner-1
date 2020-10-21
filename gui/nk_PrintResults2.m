function varargout = nk_PrintResults2(varargin)
% nk_PrintResults2 MATLAB code for nk_PrintResults2.fig
%      nk_PrintResults2, by itself, creates a new nk_PrintResults2 or raises the existing
%      singleton*.
%
%      H = nk_PrintResults2 returns the handle to a new nk_PrintResults2 or the handle to
%      the existing singleton*.
%
%      nk_PrintResults2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in nk_PrintResults2.M with the given input arguments.
%
%      nk_PrintResults2('Property','Value',...) creates a new nk_PrintResults2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nk_PrintResults2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nk_PrintResults2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nk_PrintResults2

% Last Modified by GUIDE v2.5 12-Jan-2020 17:58:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nk_PrintResults2_OpeningFcn, ...
                   'gui_OutputFcn',  @nk_PrintResults2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before nk_PrintResults2 is made visible.
function nk_PrintResults2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nk_PrintResults2 (see VARARGIN)

% Choose default command line output for nk_PrintResults2
handles.output = hObject;
handles.pnStartup.Position = handles.pnBinary.Position;
% Update handles structure
guidata(hObject, handles);
% This sets up the initial plot - only do when we are invisible
% so window can get raised using nk_PrintResults2.
if strcmp(get(hObject,'Visible'),'off')
    axes(handles.axes1); cla
    plot(rand(5));
end
set(gcf, 'WindowButtonMotionFcn', @hoverCallback)
set(gcf,'Units','Normalized','Outerposition',[0 0 0.9 0.9]);
handles.lbStartup.String = 'Initializing NM Results Viewer ...';
handles.colpt               = 'o';
handles.colptin             = getNMcolors;
handles.colpl               = '.';
handles.colpx               = '*';
handles.TitleSize           = 20;
handles.TitleWeight         = 'bold';
handles.AxisLabelSize       = 17;
handles.AxisLabelWeight     = 'bold';
handles.AxisTickSize        = 15;
handles.AxisTickWeight      = 'bold';
handles.AxisLineWidth       = 1.5;
handles.ZeroLineWidth       = 1.5;
handles.ErrorMarkerSize     = 12;
handles.ErrorMarkerWidth    = 0.25;
handles.DataMarkerSize      = 10;
handles.DataMarkerWidth     = 2.0;
handles.DataMissMarkerSize  = 12;
handles.DataMissMarkerWidth = 1;
handles.RocAxisLabelSize    = 0.15;
handles.RocAxisLabelWeight  = 'bold';
handles.RocAxisTickSize     = 0.15;
handles.RocAxisTickWeight   = 'bold';
handles.RocAxisLineWidth    = 1.5;
handles.PieAxisTickSize     = 0.15;
handles.PieAxisTickWeight   = 'bold';
handles.LegendFontSize      = 10;
handles.axes1pos_orig       = [ 0.05807522123893805 ...
                                0.08474576271186443 ...
                                0.36 ...
                                0.77 ];
handles.axes1pos_alt        = [ 0.05807522123893805 ...
                                0.39588377723970947 ...
                                0.36 ...
                                0.47 ];
set(handles.axes1,'Position', handles.axes1pos_orig);                            
set(handles.axes1, 'FontUnits','normalized','FontSize', 0.03, ...
            'FontWeight', handles.AxisTickWeight, ...
            'LineWidth', handles.AxisLineWidth);
set(handles.axes38, ...%'FontSize', handles.AxisTickSize, ...
            'FontWeight', handles.AxisTickWeight, ...
            'LineWidth', handles.AxisLineWidth);
set(handles.axes2,  'FontUnits','normalized','FontSize', 0.06, ...
            'FontWeight', handles.RocAxisTickWeight, ...
            'LineWidth', handles.RocAxisLineWidth);
set(handles.axes3, ...%'FontSize', handles.PieAxisTickSize, ...
            'FontWeight', handles.PieAxisTickWeight);
set(handles.axes4,...%'FontSize', handles.PieAxisTickSize, ...
            'FontWeight', handles.PieAxisTickWeight);
set(handles.axes17, 'FontUnits','normalized','FontSize', 0.03, ...
            'FontWeight', handles.AxisTickWeight, ...
            'LineWidth', handles.AxisLineWidth);      
        
%% Load analysis
NM = []; analind = 1; sz = get(0,'ScreenSize'); handles.screensize = sz(3:4);
if handles.screensize(1) < 800
    handles.FontSizeDynamic = 9;
elseif handles.screensize(1) < 1024
    handles.FontSizeDynamic = 10;
elseif handles.screensize(1) < 1280
    handles.FontSizeDynamic = 10.5;
elseif handles.screensize(1) < 1400
    handles.FontSizeDynamic = 11;
else
    handles.FontSizeDynamic = 11.5;
end
nVarIn = length(varargin);
for i = 1:nVarIn
   if strcmpi(varargin{i}, 'AnalysisIndex')
       analind = varargin{i+1};
   elseif strcmpi(varargin{i}, 'NM')
       NM = varargin{i+1};
   end
end

if isempty(NM), 
    try 
        NM = evalin('base', 'NM'); 
    catch
        errordlg('No NM structure found in MATLAB workspace!')
    end
end

% Read-in NM structure
handles.NM    = NM; clear NM
guidata(handles.figure1,handles)

% Populate Analysis Selector
handles = load_selAnalysis(handles);

handles.curranal = analind;

if size(handles.NM.label,2) > 1,
    handles.multilabel = true;
    handles.curlabel= get(handles.selLabel,'Value');
    handles.selLabel.String = handles.NM.labelnames;
else
    handles.multilabel = false;
    handles.curlabel = 1;
end
handles.oocvview = false;
handles.lbStartup.String = 'Retrieve OOCV results if available ...';
handles.OOCVinfo = nk_GetOOCVInfo(handles.NM,'analysis');
handles.lbStartup.String = 'Retrieve Training/CV results ...';
handles = perf_display(handles);
guidata(handles.figure1,handles);
handles.pnStartup.Visible='off';
if isfield(handles,'ExportPredictionsM')
   for i=1:numel(handles.ExportPredictionsM)
        handles.ExportPredictionsM(i).Callback = {@MenuItemPredictions_Callback, hObject};
        handles.ExportPerformanceM(i).Callback = {@MenuItemPerformance_Callback, hObject};
   end
   handles.ExportPredictionsDlg.Callback = {@MenuItemPerfTabulatur_Callback, hObject};
end

if handles.n_analyses>1
    handles.CompModels.Enable = 'on';
else
    handles.CompModels.Enable = 'off';
end

if isfield(handles,'ExportFeaturesM')
    for i=1:numel(handles.ExportFeaturesM)
        handles.ExportFeaturesM(i).Callback = {@MenuItemFeatures_Callback, hObject};
    end
end
handles.pnStartup.Visible='off';

% --- Outputs from this function are returned to the command line.
function varargout = nk_PrintResults2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

%popup_sel_index = get(handles.popupmenu1, 'Value');

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles = display_main(handles);
guidata(handles.figure1,handles)
                            
% --- Executes during object creation, after setting all properties.
function selModelMeasures_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selModelMeasures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in selYaxis.
function selYaxis_Callback(hObject, eventdata, handles)
% hObject    handle to selYaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selYaxis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selYaxis

handles = perf_display(handles);
guidata(handles.figure1, handles);

function txtPretestProb_Callback(hObject, eventdata, handles)
% hObject    handle to txtPretestProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPretestProb as text
%        str2double(get(hObject,'String')) returns contents of txtPretestProb as a double
pretest = str2double(char(get(hObject,'String')));
h_class = get(handles.popupmenu1,'Value');
display_piecharts(handles, pretest);

% --- Executes during object creation, after setting all properties.
function txtPretestProb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPretestProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in cmdUpdateX.
function cmdUpdateX_Callback(hObject, eventdata, handles)
% hObject    handle to cmdUpdateX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    % Get alternative X values from workspace
    Xvalues = evalin('base',char(get(handles.txtXAxisName,'String')));
    if ~isequal(size(Xvalues),size(handles.labels))
        errordlg('Wrong size of X value vector')
    else
        switch handles.modeflag
            case 'classification'
                for h = 1:handles.nclass
                    handles.BinClass{h}.Xaxis = Xvalues(handles.BinClass{h}.ind);
                    handles.BinClass{h}.Xaxis = handles.BinClass{h}.Xaxis(handles.BinClass{h}.sortind);
                end
            case 'regression'
                handles.Regr.Xaxis = Xvalues;
        end
        handles.XaxisName = char(get(handles.txtXAxisName,'String'));
        handles = display_main(handles);
        guidata(handles.figure1, handles);
    end
catch
    errordlg(sprintf('Could not load workspace variable %s not found!',char(get(handles.txtXAxisName,'String'))));
end

% --- Executes on button press in cmdResetX.
function cmdResetX_Callback(hObject, eventdata, handles)
% hObject    handle to cmdResetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'XaxisName')
    switch handles.modeflag
        case 'classification'
            for h = 1:handles.nclass
                handles.BinClass{h} = rmfield(handles.BinClass{h},'Xaxis');
            end
        case 'regression'
            handles.Regr = rmfield(handles.Regr,'Xaxis');
    end
    handles = rmfield(handles,'XaxisName');
    handles = display_main(handles);
    guidata(handles.figure1,handles);
end

function txtBinarize_Callback(hObject, eventdata, handles)
% hObject    handle to txtBinarize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBinarize as text
%        str2double(get(hObject,'String')) returns contents of txtBinarize as a double
handles = binarize_regr(handles);
guidata(handles.figure1,handles);

% --- Executes on selection change in selAnalysis.
function selAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to selAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selAnalysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selAnalysis

% Load analysis with pre-specified handle into main figure memory

handles = perf_display(handles);
guidata(handles.figure1,handles);

% --- Executes on selection change in selModelMeasures.
function selModelMeasures_Callback(hObject, eventdata, handles)
% hObject    handle to selModelMeasures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selModelMeasures contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selModelMeasures
handles = display_main(handles);
guidata(handles.figure1,handles);

function hoverCallback(src, evt)
    
    %handles = guidata(gcf);
    %handles=[];
    % Grab the x & y axes coordinate where the mouse is
    axesHdl = get(gcf,'CurrentAxes');
    mousePoint = get(axesHdl, 'CurrentPoint');
    figdata = get(axesHdl,'UserData');
    mousePoint2 = get(gcf,'CurrentPoint');
    if ~isempty(figdata)
        
        mouseX = mousePoint(1,1);
        mouseY = mousePoint(1,2);

        % Compare where data and current mouse point to find the data point
        % which is closest to the mouse point
        distancesToMouse = hypot(figdata.x - mouseX, figdata.y - mouseY);
        [~, indpat] = min(abs(distancesToMouse));

        % If the distance is less than some threshold, set the text
        % object's string to show the data at that point.
        xrange = nk_Range(get(axesHdl, 'Xlim'),2);
        yrange = nk_Range(get(axesHdl, 'Ylim'),2);
        if abs(mouseX - figdata.x(indpat)) < 0.035*xrange && abs(mouseY - figdata.y(indpat)) < 0.035*yrange
            figdata.textHdl.String = {figdata.patterntext{indpat}};
            figdata.hPanel.Position(1) = mousePoint2(1);
            figdata.hPanel.Position(2) = mousePoint2(2);% * (figdata.pnpos(4)/figdata.figpos(4));
            figdata.hPanel.Visible='on';
            figdata.textHdl.Visible='on';
            %uistack(figdata.hPanel,'top')
        else
            figdata.hPanel.Visible='off';
        end
        set(axesHdl,'UserData',figdata);
        
    end
    
    % --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
export_fig(handles.pnBinary,'-painters','Classification_Plot.png')
%printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Quit ' get(handles.figure1,'Name') '?'],...
                     ['Quit' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --- Executes during object creation, after setting all properties.
function txtBinarize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBinarize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selLabel.
function selLabel_Callback(hObject, eventdata, handles)
% hObject    handle to selLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selLabel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selLabel
handles = perf_display(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in selOneVsAll_Info.
function selOneVsAll_Info_Callback(hObject, eventdata, handles)
% hObject    handle to selOneVsAll_Info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selOneVsAll_Info contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selOneVsAll_Info

%% Display ROC

handles = sel_onevsone(handles, hObject);
handles = display_multiclassplot(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selOneVsAll_Info_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selOneVsAll_Info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtXAxisName_Callback(hObject, eventdata, ~)
% hObject    handle to txtXAxisName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtXAxisName as text
%        str2double(get(hObject,'String')) returns contents of txtXAxisName as a double

% --- Executes during object creation, after setting all properties.
function txtXAxisName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtXAxisName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function selAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function selYaxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selYaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cmdResetXBinary.
function cmdResetXBinary_Callback(hObject, eventdata, handles)
% hObject    handle to cmdResetXBinary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function txtXAxis_Callback(hObject, eventdata, handles)
% hObject    handle to txtXAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtXAxis as text
%        str2double(get(hObject,'String')) returns contents of txtXAxis as a double

% --- Executes during object creation, after setting all properties.
function txtXAxis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtXAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cmdUpdatePies.
function cmdUpdatePies_Callback(hObject, eventdata, ~)
% hObject    handle to cmdUpdatePies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function ExportMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ExportMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function AppearMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to AppearMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function axes17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes17


% --- Executes on selection change in selModality.
function selModality_Callback(hObject, eventdata, handles)
% hObject    handle to selModality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selModality contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selModality

handles = display_visual(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selModality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selModality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selVisMeas.
function selVisMeas_Callback(hObject, eventdata, handles)
% hObject    handle to selVisMeas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selVisMeas contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selVisMeas

handles = display_visual(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selVisMeas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selVisMeas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in selPager.
function selPager_Callback(hObject, eventdata, handles)
% hObject    handle to selPager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selPager contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selPager

handles = display_visual(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selPager_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selPager (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in tglSortFeat.
function tglSortFeat_Callback(hObject, eventdata, handles)
% hObject    handle to tglSortFeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglSortFeat

handles = display_visual(handles);
guidata(handles.figure1,handles);

% --- Executes on button press in cmdExportFeats.
function cmdExportFeats_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportFeats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

export_features(handles, 0);

% --- Executes on button press in cmdExportScores.
function cmdExportScores_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportScores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

export_scores(handles, 0);

% --- Executes on button press in cmdMetricExport.
function cmdMetricExport_Callback(hObject, eventdata, handles)
% hObject    handle to cmdMetricExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

export_performance(handles, 0);

% --- Executes on selection change in selModal.
function selModal_Callback(hObject, eventdata, handles)
% hObject    handle to selModal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selModal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selModal

handles = perf_display(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selModal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selModal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cmdSubgroupSet.
function cmdSubgroupSet_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSubgroupSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    % Get alternative X values from workspace
    I = evalin('base',char(get(handles.txtSubgroupSet,'String'))); [nI, mI] = size(I);
    if mI>1, errordlg('You entered a matrix, but a vector is required'); end
    if ~islogical(I), errordlg('A logical vector is required.'); end
    if handles.oocvview
        nC = numel(handles.OOCVinfo.Analyses{handles.curranal}.cases{handles.oocvind});
    else
        nC = size(handles.NM.label,1);
    end
    if nI ~= nC, errordlg(sprintf('The logical subgroup vector must have %g entries.',nC)); end
    handles.SubIndex = I;
    handles = perf_display(handles);
    guidata(handles.figure1,handles);
catch ERR
    errordlg(sprintf('Error: %s!',ERR.message))
end


% --- Executes on button press in cmdSubgroupClear.
function cmdSubgroupClear_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSubgroupClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'SubIndex'), handles = rmfield(handles,'SubIndex'); end
handles = perf_display(handles);
guidata(handles.figure1,handles);

function txtSubgroupSet_Callback(hObject, eventdata, handles)
% hObject    handle to txtSubgroupSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSubgroupSet as text
%        str2double(get(hObject,'String')) returns contents of txtSubgroupSet as a double


% --- Executes during object creation, after setting all properties.
function txtSubgroupSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSubgroupSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdExportAxes1.
function cmdExportAxes1_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportAxes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

titl = handles.selYaxis.String{handles.selYaxis.Value};
if strcmp(handles.axes38.Visible,'on')
    h_md_anal = [h_src handles.axes38];
else
    h_md_anal = [];
end
copyobj_subplot([handles.axes1 handles.legend_classplot h_md_anal], titl);

function [h_new, targ] = copyobj_subplot(obj, titl, targ, pos)

if ~exist('targ','var') || isempty(targ),
    targ = figure;
end

if ~exist('pos','var')
    pos = [0.12 0.15 0.85 0.7];
end

h_new = copyobj(obj,targ);
h_new(1).Position = pos; 
h_new(1).Title.String = titl;


% --- Executes on button press in cmdExportAxes2.
function cmdExportAxes2_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportAxes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
titl = 'ROC analysis';
h_ax = copyobj_subplot(handles.axes2, titl, [], [0.2 0.15 0.6 0.75] );
h_ax.Title.FontWeight = 'bold';

% --- Executes on button press in cmdExportPies.
function cmdExportPies_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportPies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prepare the pie chart descriptors
switch handles.modeflag
    case 'classification'
        h_class = get(handles.popupmenu1,'Value');
        h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
        h_classlist     = get(handles.popupmenu1,'String');
        if strcmpi(h_classlist{h_class},'Multi-group classifier') && h_onevsall_val > 1
            colind = h_onevsall_val-1;
            groupnames = {handles.NM.groupnames{colind}, 'REST'};
            cl1 = handles.colptin{colind};
            cl2 = [0.6 0.6 0.6];
        else
             groupnames = handles.BinClass{h_class}.groupnames;
            cl1 = handles.colptin(handles.BinClass{h_class}.groupind(1),:);
            cl2 = handles.colptin(handles.BinClass{h_class}.groupind(2),:);
        end
    case 'regression'
        thesh = get(handles.txtBinarize,'String');
        groupnames{1} = sprintf('Group 1 (>=%s)',thesh);
        groupnames{2} = sprintf('Group 2 (< %s)',thesh);
        cl1 = [0.5 0.5 0.5]; cl2 = [0.5 0.5 0.5];
end

hpie = copyobj_subplot(handles.axes3, [], [], [0.2 0.15 0.6 0.75] );
hpie.Title.String = sprintf('Prognostic gain: %s',groupnames{1});
hpie.Title.Position(2)=1.2;
hp      = findobj(hpie, 'Type', 'patch');
if numel(hp)==2,
    set(hp(2), 'FaceColor', cl1, 'FaceAlpha', 0.3);
end
set(hp(1), 'FaceColor', cl1, 'FaceAlpha', 1);


hpie = copyobj_subplot(handles.axes4, [], [], [0.2 0.15 0.6 0.75] );
hpie.Title.String = sprintf('Prognostic gain: %s',groupnames{2});
hpie.Title.Position(2)=1.2;
hp      = findobj(hpie, 'Type', 'patch');
if numel(hp)==2,
    set(hp(2), 'FaceColor', cl2,'FaceAlpha', 0.3);
end
set(hp(1), 'FaceColor', cl2, 'FaceAlpha', 1);

% --- Executes on button press in cmdExportCobWeb.
function cmdExportCobWeb_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportCobWeb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hcobweb = copyobj_subplot(handles.axes5, [], [], [0.2 0.10 0.6 0.75] );

% --- Executes on button press in cmdExportModelPerfVec.
function cmdExportModelPerfVec_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportModelPerfVec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
titl = 'Model performance across parameter combinations';
[h_ax, h_fig] = copyobj_subplot([handles.axes17 handles.legend_modelperf], titl, [], [0.15 0.15 0.75 0.75] );
h_ax(1).Title.FontWeight = 'bold';
h_fig.Units = 'normalized';
h_fig.Position = [0.1 0.1 0.5 0.7];

% --- Executes on button press in cmdExportModelPerfC2.
function cmdExportModelPerfC2_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportModelPerfC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
titl = 'Mean model performance across CV2 partitions';
h_ax = copyobj_subplot(handles.axes37, titl, [], [0.2 0.25 0.6 0.6] );
h_ax.Title.FontWeight = 'bold';
xlabel('Folds')
h_c=colorbar(h_ax,'location','SouthOutside');
h_c.Position(2) = 0.07; h_c.Position(4) = .05;

% --- Executes on button press in cmdExportModelPerfC1.
function cmdExportModelPerfC1_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportModelPerfC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
titl = 'Mean model performance across CV1 partitions';
h_ax = copyobj_subplot(handles.axes35, titl, [], [0.2 0.25 0.6 0.6] );
h_ax.Title.FontWeight = 'bold';
xlabel('Folds')
h_c = colorbar(h_ax,'location','SouthOutside');
h_c.Position(2) = 0.07; h_c.Position(4) = .05;


% --- Executes on button press in cmdExportAxes20.
function cmdExportAxes20_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportAxes20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
titl = 'Confusion matrix';
h_ax = copyobj_subplot(handles.axes20, titl, [], [0.2 0.15 0.6 0.75] );
h_ax.Title.FontWeight = 'bold';
colormap(flipud(gray));

% --- Executes on button press in tglSort.
function tglSort_Callback(hObject, eventdata, handles)
% hObject    handle to tglSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglSort
switch hObject.Value
    case 0
        hObject.String = 'Sort -> Values';
    case 1
        hObject.String = 'Sort -> Cases';
end
handles = perf_display(handles);
guidata(handles.figure1,handles);

% --- Executes on button press in cmdExportFeatFig.
function cmdExportFeatFig_Callback(hObject, eventdata, handles)
% hObject    handle to cmdExportFeatFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if handles.sel

titl = handles.selVisMeas.String{handles.selVisMeas.Value};
h_ax = copyobj_subplot(handles.axes33, titl, [], [0.2 0.15 0.6 0.75] );
h_ax.Title.FontWeight = 'bold';

function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selModelMeasures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in selCVoocv.
function selCVoocv_Callback(hObject, eventdata, handles)
% hObject    handle to selCVoocv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selCVoocv contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selCVoocv

handles = display_main(handles);
guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selCVoocv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selCVoocv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in selCase.
function selCase_Callback(hObject, eventdata, handles)
% hObject    handle to selCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selCase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selCase
selCase = hObject.String{hObject.Value};
AxesData = handles.axes1.UserData;
I = find(~cellfun(@isempty,strfind(AxesData.cases,selCase)));
x1 = AxesData.x(I);
y1 = AxesData.y(I);
figpos = dsxy2figxy(handles.axes1,x1,y1);
axes(handles.axes1); 
%[figx,figy] = dsxy2figxy(handles.axes1,[x1 y1],[x2 y2]);
if isfield(handles,'caseplot'), delete(handles.caseplot); end
handles.caseplot = plot(x1,y1,'ko','MarkerSize',20, 'LineWidth',1.5);

guidata(handles.figure1,handles);

% --- Executes during object creation, after setting all properties.
function selCase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function selCase_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to selCase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selCase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selCase

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over cmdExportScores.
function cmdExportScores_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to cmdExportScores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function ExportPredictionsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ExportPredictionsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function ExportPerformanceMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ExportPerformanceMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function MenuItemPerfTabulatur_Callback(hObject, eventdata, figure_handle)

handles = guidata(figure_handle);
display_tablegenerator(handles);

function MenuItemPerformance_Callback(hObject, eventdata, figure_handle)

handles = guidata(figure_handle);
PrintFromMenuClick(hObject, handles, 'export_performance')

function MenuItemPredictions_Callback(hObject, eventdata, figure_handle)

handles = guidata(figure_handle);
PrintFromMenuClick(hObject, handles, 'export_scores')

function MenuItemFeatures_Callback(hObject, eventdata, figure_handle)

handles = guidata(figure_handle);
PrintFromMenuClick(hObject, handles, 'export_features')

function PrintFromMenuClick(hObject, handles, ExportFunc)

curranal = handles.popupmenu1.Value;
if strcmp(hObject.Tag,'All')
    for i=1:numel(handles.NM.analysis)
        try
            analind = i;
            handles.curranal = analind;
            handles = switch_analysis(handles); 
            fprintf('\nProcessing analysis %g [ ID: %s ]',i, handles.NM.analysis{i}.id);
            feval(ExportFunc,handles,1);
            handles.curranal = curranal;  handles = switch_analysis(handles);
            guidata(handles.figure1,handles);
        catch
            warning('\nProcessing of analysis %g [ ID: %s ] failed',i, handles.NM.analysis{i}.id);
        end
    end
else
    
    analind = str2double(hObject.Tag);
    handles.curranal = analind;
    handles = switch_analysis(handles); 
    fprintf('\nProcessing analysis %g [ ID: %s ]',analind, handles.NM.analysis{analind}.id);
    feval(ExportFunc,handles,1);
    handles.curranal = curranal;  handles = switch_analysis(handles);
    guidata(handles.figure1,handles);
end   


% --------------------------------------------------------------------
function CompModels_Callback(hObject, eventdata, handles)
% hObject    handle to CompModels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function CompModelsCorr_Callback(hObject, eventdata, handles)
% hObject    handle to CompModelsCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CompModelsStats_Callback(hObject, eventdata, handles)
% hObject    handle to CompModelsStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

display_comparator(handles, 'stats');

% --- Executes on selection change in selSubParam.
function selSubParam_Callback(hObject, eventdata, handles)
% hObject    handle to selSubParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selSubParam contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selSubParam

handles = display_SubParam(handles);

% --- Executes during object creation, after setting all properties.
function selSubParam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selSubParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function CompModelsVisual_Callback(hObject, eventdata, handles)
% hObject    handle to CompModelsVisual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

display_comparator(handles, 'visual');
