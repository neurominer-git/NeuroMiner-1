function varargout = nk_ItemSelector(varargin)
% NK_ITEMSELECTOR MATLAB code for nk_ItemSelector.fig
%      NK_ITEMSELECTOR, by itself, creates a new NK_ITEMSELECTOR or raises the existing
%      singleton*.
%
%      H = NK_ITEMSELECTOR returns the handle to a new NK_ITEMSELECTOR or the handle to
%      the existing singleton*.
%
%      NK_ITEMSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NK_ITEMSELECTOR.M with the given input arguments.
%
%      NK_ITEMSELECTOR('Property','Value',...) creates a new NK_ITEMSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nk_ItemSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nk_ItemSelector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nk_ItemSelector

% Last Modified by GUIDE v2.5 02-Jun-2017 16:09:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nk_ItemSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @nk_ItemSelector_OutputFcn, ...
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


% --- Executes just before nk_ItemSelector is made visible.
function nk_ItemSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nk_ItemSelector (see VARARGIN)

% Choose default command line output for nk_ItemSelector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles); 
handles.M = []; handles.Cases = []; 
if(nargin > 3)
    for index = 1:2:(nargin-3),
        if nargin-3==index, break, end
        switch lower(varargin{index})
            case 'list'
                list = varargin{index+1};
            case 'matrix'
                handles.M = varargin{index+1};
                handles.YLimCases = [min(min(handles.M,[],2)) max(max(handles.M,[],2))];
                handles.YLimFeats = [min(min(handles.M)) max(max(handles.M))];
                handles.M_plot = varargin{index+1};
                set(handles.axFeats,'XTick',[]);
                set(handles.axFeats,'YTick',[]);
            case 'cases'
                handles.Cases = varargin{index+1};
            case 'selfeats'
                handles.selFeats = varargin{index+1};
            case 'selcases'
                handles.selCases = varargin{index+1};
            case 'mode'
                handles.mode = varargin{index+1};
        end
    end
end

switch handles.mode
    case 'feats'
        handles.edit2.Enable = 'off';
        handles.cmdFiltCases.Enable = 'off';
        handles.tglOpCases.Enable = 'off';
        handles.lstCasesPool.Enable = 'off';
        handles.lstCasesSelected.Enable = 'off';
        handles.cmdAddCases.Enable = 'off';
        handles.cmdRemoveCases.Enable = 'off';
        handles.cmdAddCasesAll.Enable = 'off';
        handles.cmdRemoveCasesAll.Enable = 'off';
    case 'cases'
        handles.edit1.Enable = 'off';
        handles.cmdFiltFeats.Enable = 'off';
        handles.tglOpFeats.Enable = 'off';
        handles.lstSrc.Enable = 'off';
        handles.lstDst.Enable = 'off';
        handles.cmdAdd.Enable = 'off';
        handles.cmdRemoveItems.Enable = 'off';
        handles.cmdAddAll.Enable = 'off';
        handles.cmdRemoveAll.Enable = 'off';
end

handles = load_list(handles, list, 'feats');
handles = load_list(handles, handles.Cases, 'cases');
handles = display_matrix(handles);
handles = display_eval(handles);

figdata = DefineMatrixTextBox(handles,handles.Cases,handles.Items);
set(handles.axMat,'UserData',figdata);
set(gcf, 'WindowButtonMotionFcn', @hoverCallback);

% Update handles structure
guidata(handles.figure1, handles);

% Make the GUI modal
%set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes nk_OOCVDataIO wait for user response (see UIRESUME)
uiwait(handles.figure1);

% UIWAIT makes nk_ItemSelector wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function hoverCallback(src, evt)
    
    %handles = guidata(gcf);
    %handles=[];
    % Grab the x & y axes coordinate where the mouse is
    axesHdl = get(gcf,'CurrentAxes');
    mousePoint = get(axesHdl, 'CurrentPoint');
    figdata = get(axesHdl,'UserData');
    if ~isempty(figdata)
        
        mouseX = mousePoint(1,1);
        mouseY = mousePoint(1,2);

        % Compare where data and current mouse point to find the data point
        % which is closest to the mouse point
        distancesToMouse = hypot(figdata.x - mouseX, figdata.y - mouseY);
        [~, indpat] = min(abs(distancesToMouse(:)));

        % If the distance is less than some threshold, set the text
        % object's string to show the data at that point.
        xrange = nk_Range(get(axesHdl, 'Xlim'),2);
        yrange = nk_Range(get(axesHdl, 'Ylim'),2);
        if abs(mouseX - figdata.x(indpat)) < 0.02*xrange && abs(mouseY - figdata.y(indpat)) < 0.02*yrange
            figdata.lbCaseFeat.String =  figdata.patterntext(indpat);
            figdata.lbCaseFeat.Visible='on';
        else
            figdata.lbCaseFeat.Visible='off';
        end
        set(axesHdl,'UserData',figdata);
        
    end

function figdata = DefineMatrixTextBox(handles,cases,list)
    
nC = numel(cases);
nL = numel(list);

% Define textbox info data 
pss = cell(nC,nL);

for j = 1:nC
    for i = 1:nL
        pss{j,i} = sprintf('Case [%g]: %s  |  Feature [%g]: %s', j, cases{j}, i, list{i});
    end
end

figdata.x = repmat(1:numel(list),nC,1);
figdata.y = repmat((1:numel(cases))',1,nL);
figdata.patterntext = pss;
figdata.lbCaseFeat = handles.lbCaseFeat;

function handles = load_list(handles, list, act)

switch act
    case 'feats'
        handles.Items = list;
        if isfield(handles,'selFeats') && ~isempty(handles.selFeats) 
            selFeats = handles.selFeats; 
        else
            selFeats = true(numel(list),1); 
        end
        set(handles.lstDst, 'String', list(selFeats));
        set(handles.lstDst, 'Max', numel(list(selFeats)),'Min',0);
        if ~any(selFeats==0)
            set(handles.lstSrc, 'String', {''});
            selNFeats = [];
        else
            set(handles.lstSrc, 'String', list(~selFeats));
            selNFeats = ~selFeats;
        end
        handles.DstIndex = find(selFeats);
        handles.SrcIndex = find(selNFeats);
        
    case 'cases'
        handles.Cases = list;
        if isfield(handles,'selCases') && ~isempty(handles.selCases) 
            selCases = handles.selCases; 
        else
            selCases = true(numel(list),1); 
        end
        set(handles.lstCasesSelected, 'String', list(selCases));
        set(handles.lstCasesSelected, 'Max', numel(list(selCases)),'Min',0);
        if ~any(selCases==0)
            set(handles.lstCasesPool, 'String', {''});
            selNCases = [];
        else
            set(handles.lstCasesPool, 'String', list(~selCases));
            selNCases = ~selCases;
        end
        handles.SelectedCasesIndex = find(selCases);
        handles.AllCasesIndex = find(selNCases);
end

function handles = display_matrix(handles)

if ~isempty(handles.M_plot)
   handles.M_plot = filter_matrix(handles.M, handles);
   axes(handles.axMat); cla
   handles.mat_plot = imagesc(handles.M_plot, 'Parent', handles.axMat);
   set(handles.mat_plot,'ButtonDownFcn', {@axMat_ButtonDownFcn,handles});
   %set(handles.mat_plot,'HitTest','off');
   if isfield(handles,'selCases') && ~isempty(handles.selCases), Cases = handles.Cases(handles.selCases); else Cases = handles.Cases; end
   if isfield(handles,'selFeats') &&~isempty(handles.selFeats), Feats = handles.Items(handles.selFeats); else Feats = handles.Items; end
   figdata = DefineMatrixTextBox(handles,Cases,Feats);
   set(handles.axMat,'UserData',figdata);
end

% --- Outputs from this function are returned to the command line.
function varargout = nk_ItemSelector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure);

if isvalid(hObject)
    if isfield(handles,'selFeats'), selFeats = handles.selFeats; else, selFeats = []; end
    varargout{1} = selFeats;
    if isfield(handles,'selCases'), selCases = handles.selCases; else, selCases = []; end
    varargout{2} = selCases;
    % The figure can be deleted now
    delete(handles.figure1);
else
    varargout{1} = [];
    varargout{2} = [];
end

% --- Executes on selection change in lstDst.
function lstDst_Callback(hObject, eventdata, handles)
% hObject    handle to lstDst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstDst contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstDst


% --- Executes during object creation, after setting all properties.
function lstDst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstDst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in lstSrc.
function lstSrc_Callback(hObject, eventdata, handles)
% hObject    handle to lstSrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstSrc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstSrc


% --- Executes during object creation, after setting all properties.
function lstSrc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstSrc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ModItems(cmd, handles, act)

switch cmd
    case {'rem','remall'}
        switch act
            case 'feats'
                AddInd  = 'SrcIndex'; 
                RemInd  = 'DstIndex'; 
                LstSrc  = 'lstDst';
                LstDst  = 'lstSrc';
                T='Items';                
                TT='selFeats';
                D = 2;
            case 'cases'
                AddInd  = 'AllCasesIndex'; 
                RemInd  = 'SelectedCasesIndex'; 
                LstSrc  = 'lstCasesSelected';
                LstDst  = 'lstCasesPool';
                T='Cases';                
                TT='selCases';
                D = 1;
        end
        
    case {'add','addall'}
        switch act
            case 'feats'
                AddInd  = 'DstIndex'; 
                RemInd  = 'SrcIndex'; 
                LstSrc  = 'lstSrc'; 
                LstDst  = 'lstDst';   
                T='Items';                
                TT='selFeats';
                D = 2;
            case 'cases'
                RemInd  = 'AllCasesIndex'; 
                AddInd  = 'SelectedCasesIndex'; 
                LstDst  = 'lstCasesSelected';
                LstSrc  = 'lstCasesPool';
                T='Cases';                
                TT='selCases';
                D = 1;
        end
end

if strcmp(handles.(LstSrc).String,''); return; end

if strfind(char(cmd),'all')
    handles.(LstSrc).Value=1:numel(handles.(LstSrc).String);
end

selItems = get(handles.(LstSrc),'Value');

if ~isempty(selItems)
    switch act
        case 'feats'
            handles.(AddInd) = sort([ handles.(AddInd); handles.(RemInd)(selItems) ]);
        case 'cases'
            handles.(AddInd) = sort([ handles.(AddInd); handles.(RemInd)(selItems) ]);
    end
else
    handles.(AddInd) = sort( handles.(AddInd) );
end
handles.(RemInd)(selItems)=[];
if ~numel(handles.(AddInd))
    set(handles.(LstDst),'String',{''});
    set(handles.(LstDst),'Max',1)
    set(handles.(LstDst),'Value',1);
else
    handles.(LstDst).String = handles.(T)(handles.(AddInd));
    handles.(LstDst).Max = numel(handles.(T)(handles.(AddInd)));
end
if ~numel(handles.(RemInd))
    set(handles.(LstSrc),'String',{''});
    set(handles.(LstSrc),'Max',1);
    set(handles.(LstSrc),'Value',1);
else
    handles.(LstSrc).String = handles.(T)(handles.(RemInd));
    handles.(LstSrc).Max = numel(handles.(T)(handles.(RemInd)));
end
nS = numel(handles.(LstSrc).String);
nD = numel(handles.(LstDst).String);
if nD > 1, 
    handles.(LstDst).Value = numel(handles.(AddInd))-1 ;
else
    handles.(LstDst).Value = 1;
end    
if nS > 1, 
    handles.(LstSrc).Value = numel(handles.(RemInd))-1 ;
else
    handles.(LstSrc).Value = 1;
end  

switch act
    case 'feats'
        handles.(TT) = false(1,size(handles.M, D));
        handles.(TT)(handles.DstIndex) = true;
    case 'cases'
        handles.(TT) = false(size(handles.M, D),1);
        handles.(TT)(handles.SelectedCasesIndex) = true;
end
if ~any(handles.(TT))
    cla(handles.axMat);
    cla(handles.axFeats);
    cla(handles.axCases);
else
    handles = display_matrix(handles);
    handles = display_eval(handles);
end

% Update handles structure
guidata(handles.figure1, handles);

% --- Executes on button press in cmdOK.
function cmdOK_Callback(hObject, eventdata, handles)
% hObject    handle to cmdOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on selection change in selMeas.
function selMeas_Callback(hObject, eventdata, handles)
% hObject    handle to selMeas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selMeas contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selMeas

handles = display_matrix(handles);
handles = display_eval(handles);

% Update handles structure
guidata(handles.figure1, handles)

% --- Executes during object creation, after setting all properties.
function selMeas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selMeas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = display_eval(handles)

handles.M_plot = filter_matrix(handles.M, handles);
str = handles.selMeas.String{handles.selMeas.Value};
hold(handles.axFeats,'off');
hold(handles.axCases,'off');

switch str
    case 'mean'
        mF = nm_nanmean(handles.M_plot,1);
        mC = nm_nanmean(handles.M_plot,2);
    case 'median'
        mF = nm_nanmedian(handles.M_plot,1);
        mC = nm_nanmedian(handles.M_plot,2);
    case 'variance'
        mF = nm_nanvar(handles.M_plot,1);
        mC = nm_nanvar(handles.M_plot,2);
    case 'interquartile range'
        mF = iqr(handles.M_plot,1);
        mC = iqr(handles.M_plot,2);
    case 'identical'
        mC = zeros(size(handles.M_plot,1),1);
        mF = zeros(size(handles.M_plot,2),1);
        for i=1:size(handles.M_plot,2)
            mF(i) = numel(unique(handles.M_plot(:,i)))==1;
        end
        for i=1:size(handles.M_plot,1)
            mC(i) = numel(unique(handles.M_plot(i,:)))==1;
        end  
    case '%(nan)'
        mF = sum(isnan(handles.M_plot),1)*100/size(handles.M_plot,1);
        mC = sum(isnan(handles.M_plot),2)*100/size(handles.M_plot,2);
    case 'skewness'
        mF = skewness(handles.M_plot,1,1);
        mC = skewness(handles.M_plot,1,2);
    case 'kurtosis'
        mF = kurtosis(handles.M_plot,1,1);
        mC = kurtosis(handles.M_plot,1,2);
    case 'cross-correlation'
        IN.X = handles.M_plot;
        X = nk_PerfImputeObj(handles.M_plot,IN);
        nF = size(handles.M_plot,2); mF = zeros(1,nF);
        for i=1:nF
            ind = true(1,nF); ind(i)=false; 
            Fi = X(:,i); Fni = X(:,ind);
            mF(i) = mean(abs(nk_CorrMat(Fni,Fi)));
        end
        nC = size(handles.M_plot,1); mC = zeros(1,nC);
        for i=1:nC
            ind = true(1,nC); ind(i)=false; 
            Ci = X(i,:); Cni = X(ind,:);
            mC(i) = mean(abs(nk_CorrMat(Cni',Ci')));
        end
end

handles.feats_plot = plot(handles.axFeats,mF,'ko-','LineWidth',1,'MarkerSize', 3); 
if size(handles.M_plot,2)>1,xlim(handles.axFeats,[1 size(handles.M_plot,2)]);end
handles.cases_plot = plot(handles.axCases,mC,'ko-','LineWidth',1,'MarkerSize', 3);
if size(handles.M_plot,1)>1,xlim(handles.axCases,[1 size(handles.M_plot,1)]);end
handles.axCases.XDir = 'reverse';
hold(handles.axFeats,'on');
hold(handles.axCases,'on');

if isfield(handles,'thr_feats')
    if isfield(handles,'thrCases'), delete(handles.thrCases);end
    handles.thrFeats = plotpatch(handles.M_plot, handles.thr_feats, handles.axFeats, handles.tglOpFeats);
end
if isfield(handles,'thr_cases')
    if isfield(handles,'thrCases'), delete(handles.thrCases);end
    handles.thrCases = plotpatch(handles.M_plot, handles.thr_cases, handles.axCases, handles.tglOpCases);
end
handles.mF = mF;
handles.mC = mC;
view(handles.axCases,[-90,90]) 
handles.axCases.YDir = 'reverse';
handles.axCases.XTick = [];
handles.axFeats.XTick = [];

% Update handles structure
guidata(handles.figure1, handles)

function handle = plotpatch(M, thr, ax, tgl)

if strcmp(tgl.String,'>=')
    st = thr; hg = ax.YLim(end);
else
    st = ax.YLim(1); hg = thr;
end
handle = patch(ax, [1 size(M,1) size(M,1) 1], [ st st hg hg ], 'red', 'FaceAlpha', 0.5, 'LineStyle', 'none');

% --- Executes on button press in tglScale.
function tglScale_Callback(hObject, eventdata, handles)
% hObject    handle to tglScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglScale

if get(hObject,'Value') 
   handles.M_unscaled = handles.M;
   handles.M = nk_PerfScaleObj(handles.M,[]);
else
   handles.M = handles.M_unscaled;
end

handles = display_matrix(handles);
handles = display_eval(handles);

% Update handles structure
guidata(handles.figure1, handles)

function M = filter_matrix(M, handles)

if isfield(handles,'selFeats') && ~isempty(handles.selFeats)
    M = M(:,handles.selFeats);
end
if isfield(handles,'selCases') && ~isempty(handles.selCases)
    M = M(handles.selCases,:);
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.thr_feats = str2double(char(hObject.String));
if isfield(handles,'thrFeats'), delete(handles.thrFeats); end
set(hObject,'BackgroundColor',orange);
handles = display_eval(handles);

% Update handles structure
guidata(handles.figure1, handles)
 
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles.thr_cases = str2double(char(hObject.String));
if isfield(handles,'thrCases'), delete(handles.thrCases); end
set(hObject,'BackgroundColor',orange);
handles = display_eval(handles);
% Update handles structure
guidata(handles.figure1, handles)

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [sel_ind, thr] = eval_filter(thr, ax, lims, vals, sel_prev, op)

sel_ind = true(lims(2)-lims(2)+1,1);

switch op
    case '>='
        op = @le;
    case '<='
        op = @ge;
end

hold(ax,'on');
try 
   sel_ind = op(vals,thr); 
   if ~isempty(sel_prev), 
       sel_prev_f = find(sel_prev);
       sel_new_f = sel_prev_f(sel_ind);
       sel_ind = false(lims(2)-lims(2)+1,1);
       sel_ind(sel_new_f)=true;
   end
   if ~any(sel_ind)
       errordlg('Your filtering was too aggressive. Returning to previous setting')
       sel_ind = true(lims(2)-lims(2)+1,1);
       thr=max(vals);
   end
catch
   errordlg('Invalid input into field!','Error')
end

% --- Executes on button press in cmdFiltFeats.
function cmdFiltFeats_Callback(hObject, eventdata, handles)
% hObject    handle to cmdFiltFeats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'selFeats'), selFeats = handles.selFeats; else, selFeats = []; end
if isfield(handles,'thr_feats'), 
    thr_feats=handles.thr_feats;
else
    thr_feats = str2double(char(handles.edit1.String));
end
[handles.selFeats, handles.thr_feats] = eval_filter(thr_feats, handles.axFeats, [1 size(handles.M,2)], handles.mF, selFeats, handles.tglOpFeats.String);
set(handles.edit1,'BackgroundColor','red');

handles = display_matrix(handles);
handles = display_eval(handles);
handles = load_list(handles, handles.Items,'feats');
handles = load_list(handles, handles.Cases,'cases');
% Update handles structure
guidata(handles.figure1, handles)

% --- Executes on button press in cmdFiltCases.
function cmdFiltCases_Callback(hObject, eventdata, handles)
% hObject    handle to cmdFiltCases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'thrCases'), delete(handles.thrCases); end
if isfield(handles,'thr_cases'), 
    thr_cases=handles.thr_cases;
else
    thr_cases = str2double(char(handles.edit2.String));
end
if isfield(handles,'selCases'), selCases = handles.selCases; else, selCases = []; end
[handles.selCases, handles.thr_cases] = eval_filter(thr_cases, handles.axCases, [1 size(handles.M,1)], handles.mC, selCases, handles.tglOpCases.String);
set(handles.edit2,'BackgroundColor','red');

handles = display_matrix(handles);
handles = display_eval(handles);
handles = load_list(handles, handles.Items,'feats');
handles = load_list(handles, handles.Cases,'cases');

% Update handles structure
guidata(handles.figure1, handles)

% --- Executes on button press in cmdResetFilt.
function cmdResetFilt_Callback(hObject, eventdata, handles)
% hObject    handle to cmdResetFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'selFeats'),
    handles = rmfield(handles,'selFeats'); 
end
if isfield(handles,'selCases'),
    handles = rmfield(handles,'selCases'); 
end
if isfield(handles,'thrFeats')
    delete(handles.thrFeats); handles = rmfield(handles,'thrFeats'); handles = rmfield(handles,'thr_feats');
end
if isfield(handles,'thrCases')
    delete(handles.thrCases); handles = rmfield(handles,'thrCases'); handles = rmfield(handles,'thr_cases');
end

set(handles.edit1,'BackgroundColor','white');
set(handles.edit2,'BackgroundColor','white');
handles.edit1.String = {'0'};
handles.edit2.String = {'0'};

handles = display_matrix(handles);
handles = display_eval(handles);
handles = load_list(handles, handles.Items,'feats');
handles = load_list(handles, handles.Cases,'cases');

% Update handles structure
guidata(handles.figure1, handles)

% --- Executes on button press in cmdAdd.
function cmdAdd_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('add',handles, 'feats')

% --- Executes on button press in cmdRemoveItems.
function cmdRemoveItems_Callback(hObject, eventdata, handles)
% hObject    handle to cmdRemoveItems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('rem',handles, 'feats')

% --- Executes on button press in cmdRemoveAll.
function cmdRemoveAll_Callback(hObject, eventdata, handles)
% hObject    handle to cmdRemoveAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('remall',handles,'feats')

% --- Executes on button press in cmdAddAll.
function cmdAddAll_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAddAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('addall',handles,'feats')

function [C] = orange()
C = [1 .5 0];

% --- Executes on mouse press over axes background.
function axMat_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figdata = get(handles.axMat,'UserData');
mousePoint = get(handles.axMat, 'CurrentPoint');
mouseX = mousePoint(1,1); mouseY = mousePoint(1,2);
distancesToMouse = hypot(figdata.x - mouseX, figdata.y - mouseY);
[~, indpat] = min(abs(distancesToMouse(:)));
xrange = nk_Range(get(handles.axMat, 'Xlim'),2);
yrange = nk_Range(get(handles.axMat, 'Ylim'),2);
if abs(mouseX - figdata.x(indpat)) < 0.02*xrange && abs(mouseY - figdata.y(indpat)) < 0.02*yrange
  handles.lstDst.Value = figdata.x(indpat);
end

% --- Executes on button press in tglOpFeats.
function tglOpFeats_Callback(hObject, eventdata, handles)
% hObject    handle to tglOpFeats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglOpFeats

if get(hObject,'Value') 
   hObject.String = '<=';
else
   hObject.String = '>=';
end
if isfield(handles,'thr_feats')
    if isfield(handles,'thrFeats'), delete(handles.thrFeats);end
    handles.thrFeats = plotpatch(handles.M_plot, handles.thr_feats, handles.axFeats, handles.tglOpFeats);
    % Update handles structure
    guidata(handles.figure1, handles)
end

% --- Executes on button press in tglOpFeats.
function tglOpCases_Callback(hObject, eventdata, handles)
% hObject    handle to tglOpFeats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglOpFeats

if get(hObject,'Value') 
   hObject.String = '<=';
else
   hObject.String = '>=';
end
if isfield(handles,'thr_cases')
    if isfield(handles,'thrCases'), delete(handles.thrCases);end
    handles.thrCases = plotpatch(handles.M_plot, handles.thr_cases, handles.axCases, handles.tglOpCases);
    % Update handles structure
    guidata(handles.figure1, handles)
end


% --- Executes on button press in cmdAddCases.
function cmdAddCases_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAddCases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('add',handles,'cases')

% --- Executes on button press in cmdRemoveCases.
function cmdRemoveCases_Callback(hObject, eventdata, handles)
% hObject    handle to cmdRemoveCases (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('rem',handles,'cases')

% --- Executes on button press in cmdRemoveCasesAll.
function cmdRemoveCasesAll_Callback(hObject, eventdata, handles)
% hObject    handle to cmdRemoveCasesAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('remall',handles,'cases')

% --- Executes on button press in cmdAddCasesAll.
function cmdAddCasesAll_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAddCasesAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ModItems('addall',handles,'cases')

% --- Executes on selection change in lstCasesPool.
function lstCasesPool_Callback(hObject, eventdata, handles)
% hObject    handle to lstCasesPool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstCasesPool contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstCasesPool

% --- Executes during object creation, after setting all properties.
function lstCasesPool_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstCasesPool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lstCasesSelected.
function lstCasesSelected_Callback(hObject, eventdata, handles)
% hObject    handle to lstCasesSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstCasesSelected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstCasesSelected


% --- Executes during object creation, after setting all properties.
function lstCasesSelected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstCasesSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
