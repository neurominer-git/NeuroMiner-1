function varargout = nk_main(varargin)
% NK_MAIN MATLAB code for nk_main.fig
%      NK_MAIN, by itself, creates a new NK_MAIN or raises the existing
%      singleton*.
%
%      H = NK_MAIN returns the handle to a new NK_MAIN or the handle to
%      the existing singleton*.
%
%      NK_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NK_MAIN.M with the given input arguments.
%
%      NK_MAIN('Property','Value',...) creates a new NK_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nk_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nk_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nk_main

% Last Modified by GUIDE v2.5 04-Dec-2014 15:52:39

% Begin initialization code - DO NOT EDIT 

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nk_main_OpeningFcn, ...
                   'gui_OutputFcn',  @nk_main_OutputFcn, ...
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
end

% --- Executes just before nk_main is made visible.
function nk_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nk_main (see VARARGIN)
global NM
try
    NM = evalin('base','NM'); 
catch
    loadflag = questdlg(['No NeuroMiner structure found in workspace. ' ...
        'Load from disk ?'],'NM Startup','Yes');
    if strcmp(loadflag,'Yes'), NM = loadmat; end
end
 
% Choose default command line output for nk_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nk_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end

% --- Outputs from this function are returned to the command line.
function varargout = nk_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --------------------------------------------------------------------
function file_main_Callback(hObject, eventdata, handles)
% hObject    handle to file_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function file_load_data_Callback(hObject, eventdata, handles)
% hObject    handle to file_load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --------------------------------------------------------------------
function file_load_NM_Callback(hObject, eventdata, handles)
% hObject    handle to file_load_NM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NM
NM = loadmat;
end

% --------------------------------------------------------------------
function file_save_NM_Callback(hObject, eventdata, handles)
% hObject    handle to file_save_NM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NM
savemat(NM)
end

% --------------------------------------------------------------------
function file_quit_NM_Callback(hObject, eventdata, handles)
% hObject    handle to file_quit_NM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QuitNeuroMiner
delete(handles.figure1)
end

% --------------------------------------------------------------------
function file_cd_pwd_Callback(hObject, eventdata, handles)
% hObject    handle to file_cd_pwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
directoryname = uigetdir(pwd, 'Pick a new working directory');
if directoryname, cd(directoryname); end
end

% --------------------------------------------------------------------
function file_del_NM_Callback(hObject, eventdata, handles)
% hObject    handle to file_del_NM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearflag = questdlg('Are you sure you want delete the NM structure ?','Delete NM structure','No');
if strcmp(clearflag,'yes'), NM=[]; end

end

function NM = loadmat

[filename, pathname] = uigetfile('*.mat', 'Pick a NeuroMiner structure file');

if ~isequal(filename,0) && ~isequal(pathname,0)
    fprintf('\nLoading %s ... please wait',filename)
    load(fullfile(pathname,filename))
    if ~exist('NM','var')
        retryflag = ...
            nk_input([filename ' does not contain a valid NeuroMiner structure! Retry ?'], ...
            0, 'yes|no',[1,0],1);
        if retryflag, loadmat, end
    else
        cd(pathname)
    end
else
    NM = [];
end

end

function savemat(NM)

[filename, pathname] = uiputfile('*.mat', 'Save NeuroMiner structure file');

if ~isequal(filename,0) && ~isequal(pathname,0)
    fprintf('\nSaving %s ... please wait',filename)
    save(fullfile(pathname,filename),'NM')
end

end

function QuitNeuroMiner
global NM

fprintf('\nGood bye.\n')
assignin('base', 'NM', NM)
clearvars -global NM

end