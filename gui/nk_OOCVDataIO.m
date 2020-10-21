function varargout = nk_OOCVDataIO(varargin)
% NK_OOCVDATAIO MATLAB code for nk_OOCVDataIO.fig
%      NK_OOCVDATAIO by itself, creates a new NK_OOCVDATAIO or raises the
%      existing singleton*.
%
%      H = NK_OOCVDATAIO returns the handle to a new NK_OOCVDATAIO or the handle to
%      the existing singleton*.
%
%      NK_OOCVDATAIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NK_OOCVDATAIO.M with the given input arguments.
%
%      NK_OOCVDATAIO('Property','Value',...) creates a new NK_OOCVDATAIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nk_OOCVDataIO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nk_OOCVDataIO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nk_OOCVDataIO

% Last Modified by GUIDE v2.5 27-Jan-2017 22:50:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nk_OOCVDataIO_OpeningFcn, ...
                   'gui_OutputFcn',  @nk_OOCVDataIO_OutputFcn, ...
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

% --- Executes just before nk_OOCVDataIO is made visible.
function nk_OOCVDataIO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nk_OOCVDataIO (see VARARGIN)

% Choose default command line output for nk_OOCVDataIO
handles.output = [];

% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
             case 'title'
                set(hObject, 'Name', varargin{index+1});
             case 'list'
                 if ~isempty(varargin{index+1})
                    for i=1:numel(varargin{index+1})
                        handles.Data.Items{i} = varargin{index+1}{i};
                    end
                 else
                    handles.Data.Items=[];
                 end
                 guidata(handles.figure1,handles);
                 print_dataitems(handles)
              case 'mode'
                 switch varargin{index+1}
                     case 'select'
                         bt_str = 'Select dataset'; vis_str = 'off';
                         set(handles.lstData,'Position',[0.034 0.05 0.93 0.82]);
                         handles.mode = 0;
                     case 'modify'
                         bt_str = 'Add/Modify Data in Dataset'; vis_str = 'on';
                         set(handles.lstData,'Position',[0.034 0.17 0.93 0.70]);
                         handles.mode = 1;
                 end
                 set(handles.cmdSelect,'String',bt_str);
                 set(handles.cmdDelete,'Visible',vis_str);
                 set(handles.txtNewData,'Visible',vis_str);
                 set(handles.cmdCreateNew,'Visible',vis_str);
                 set(handles.chkLabelKnown,'Visible',vis_str);
                 set(handles.chkCalibrationDataAvail,'Visible',vis_str);
                 set(handles.cmdSave,'Visible',vis_str);
        end
    end
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes nk_OOCVDataIO wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = nk_OOCVDataIO_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in cmdSelect.
function cmdSelect_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Data.SelItem = get(handles.lstData,'Value');
handles.output = handles.Data;

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = [];

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    


% --- Executes on selection change in lstData.
function lstData_Callback(hObject, eventdata, handles)
% hObject    handle to lstData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstData contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstData


% --- Executes during object creation, after setting all properties.
function lstData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function txtNewData_Callback(hObject, eventdata, handles)
% hObject    handle to txtNewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNewData as text
%        str2double(get(hObject,'String')) returns contents of txtNewData as a double


% --- Executes during object creation, after setting all properties.
function txtNewData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdCreateNew.
function cmdCreateNew_Callback(hObject, eventdata, handles)
% hObject    handle to cmdCreateNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.Data,'Items'), l = numel(handles.Data.Items)+1; else l=1; end
handles.Data.Items{l} = struct('desc', get(handles.txtNewData,'String'), ...
                                'date',datestr(now), ...
                                'n_subjects_all',0, ...
                                'labels_known', get(handles.chkLabelKnown,'Value'), ...
                                'os', sprintf('%s', system_dependent('getos')));
if isfield(handles.Data,'NewItemIndex'), lx = numel(handles.Data.NewItemIndex)+1; else lx = 1; end
handles.Data.NewItemIndex(lx) = l;
print_dataitems(handles)
guidata(handles.figure1,handles);

% --- Executes on button press in cmdDelete.
function cmdDelete_Callback(hObject, eventdata, handles)
% hObject    handle to cmdDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selItem = get(handles.lstData,'Value');
lstItems = get(handles.lstData,'String');
delfl = questdlg(sprintf('Are you sure you want to delete dataset:\n%s',lstItems{selItem}));
if strcmp(delfl,'Yes')
%     handles.Data.Items(selItem) = [];
%     if isfield(handles.Data,'NewItems')
%         ind = find(handles.Data.NewItems, selItem);
%         handles.Data.NewItems(ind)=[];
%     end

    handles.output.delete = selItem;
    guidata(handles.figure1,handles);
    uiresume(handles.figure1)

%     % Update handles structure
%     guidata(hObject, handles);
%     
%     print_dataitems(handles,selItem-1);
end

function toggle_controls(handles)

if ~isfield(handles,'Data') || ~isfield(handles.Data,'Items') || isempty(handles.Data.Items);
    tglstr = 'off';
else
    tglstr = 'on';
end

set(handles.cmdSelect,'Enable',tglstr);
set(handles.cmdSave,'Enable',tglstr);
set(handles.cmdDelete,'Enable',tglstr);
set(handles.lstData,'Enable',tglstr);

function print_dataitems(handles, value)
if ~exist('value','var') || isempty(value), value = numel(handles.Data.Items); end
if ~isempty(handles.Data.Items)
    for i=1:numel(handles.Data.Items)
        if isfield(handles.Data.Items{i},'labels_known')
            if handles.Data.Items{i}.labels_known
                lbknown_str = '; labels known: yes';
            else
                lbknown_str = '; labels known: no';
            end
        else
            lbknown_str = '';
        end
        if ~isfield(handles.Data.Items{i},'n_subjects_all')
            n_subjects_all=0;
        else
            n_subjects_all= handles.Data.Items{i}.n_subjects_all;
        end
        datastr{i} = sprintf('[%3i] %s (%g cases); date: %s%s',i,handles.Data.Items{i}.desc,n_subjects_all,handles.Data.Items{i}.date, lbknown_str);
    end
    set(handles.lstData, 'String', datastr);
    handles.lstData.Value=value;
else
    set(handles.lstData, 'String', 'NO DATA AVALABLE');
    handles.lstData.Value=1;
end
 toggle_controls(handles)

% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selItem = get(handles.lstData,'Value');
%lstItems = get(handles.lstData,'String');
OOCV = handles.Data.Items{selItem};
filename = 'NM_OOCV_';
uisave('OOCV',filename);

% --- Executes on button press in chkLabelKnown.
function chkLabelKnown_Callback(hObject, eventdata, handles)
% hObject    handle to chkLabelKnown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkLabelKnown


% --- Executes on button press in chkCalibrationDataAvail.
function chkCalibrationDataAvail_Callback(hObject, eventdata, handles)
% hObject    handle to chkCalibrationDataAvail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkCalibrationDataAvail


% --------------------------------------------------------------------
function uiLoadData_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uiLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uiNewData_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uiNewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdLoad.
function cmdLoad_Callback(hObject, eventdata, handles)
% hObject    handle to cmdLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n = numel(handles.Data.Items);

[filename,pathname] = uigetfile('*.mat','Load Independent Test Data','MultiSelect','on');
if iscell(filename)
    nF = numel(filename);
else
    nF=1;
end
if filename 
    for i=1:nF
        if iscell(filename)
            pth = fullfile(pathname,filename{i});
            load(pth)
        else
            pth = fullfile(pathname,filename);
            load(pth)
        end
        handles.Data.Items{n+i} = OOCV;
    end
end
guidata(handles.figure1,handles);
print_dataitems(handles)
