function varargout = nk_results(varargin)
% NK_RESULTS MATLAB code for nk_results.fig
%      NK_RESULTS, by itself, creates a new NK_RESULTS or raises the existing
%      singleton*.
%
%      H = NK_RESULTS returns the handle to a new NK_RESULTS or the handle to
%      the existing singleton*.
%
%      NK_RESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NK_RESULTS.M with the given input arguments.
%
%      NK_RESULTS('Property','Value',...) creates a new NK_RESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nk_results_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nk_results_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nk_results

% Last Modified by GUIDE v2.5 15-Dec-2014 16:44:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nk_results_OpeningFcn, ...
                   'gui_OutputFcn',  @nk_results_OutputFcn, ...
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

% --- Executes just before nk_results is made visible.
function nk_results_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nk_results (see VARARGIN)

% Choose default command line output for nk_results
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using nk_results.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes nk_results wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nk_results_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        print_classification(NM, analind, h, GraphType)
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end

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

modeflag = 'classification';

switch modeflag

    case 'classification'
        lstmnu = {'Classification plot', ... 
                    'ROC plot', ...
                    'CV1 / CV2 performance plot', ...
                    'Model complexity plot', ...
                    'Subject prediction list'};
                
    case 'regression'
        lstmnu = {'Regression plot', ...
                     'CV1 / CV2 performance plot', ...
                    'Model complexity plot', ...
                    'Subject prediction list'};

end

set(hObject, 'String', lstmnu);


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

strmnu = {'Means without error bars', ...
        'Means with standard deviation', ...
        'Means with 95% confidence intervals', ...
        'Probabilities'};

set(hObject, 'String', strmnu);

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function print_classification(NM, analind, h, GraphType)

SVM = NM.analysis{analind}.params.SVM;
CV = NM.analysis{analind}.params.SVM;
nclass = numel(NM.analysis{analind}.BinClass);
BinClass = NM.analysis{analind}.Binclass{h};
class = CV.class{1,1}{h};
cases = NM.cases;
label = NM.label; [label,s] = sort(label,'ascend');

colpt   = {'bo','ro','go','mo','co','ko'};
colptin = {'b','r','g','m','c','k'};

%% Formatting
TitleSize           = 20;
TitleWeight         = 'bold';
AxisLabelSize       = 18;
AxisLabelWeight     = 'bold';
AxisTickSize        = 18;
AxisTickWeight      = 'bold';
AxisLineWidth       = 1.5;
ZeroLineWidth       = 1.5;
ErrorMarkerSize     = 12;
ErrorMarkerWidth    = 1;
DataMarkerSize      = 12;
DataMarkerWidth     = 2.5;
DataMissMarkerSize  = 18;
DataMissMarkerWidth = 3;

groupind = class.groups;
   
% Setup some variables for convenience
groupnames = regexp(class.groupdesc,' vs ','split');

ind     = BinClass.index_predictions(s);

ha = gca; set(ha,'FontSize',AxisTickSize,'FontWeight', AxisTickWeight,'LineWidth', AxisLineWidth); 
box on; hold on

labelh = zeros(size(label,1),1);

if numel(groupind) == 2
    ind1 = label == groupind(1); ind2 = label == groupind(2);
    labelh(ind1,1) = 1; labelh(ind2,1) = -1;
else
    ind1 = label == groupind(1); 
    labelh(ind1,1) = 1; labelh(~ind,1) = -1;
end

if find(groupind==0), groupind = groupind + 1; end
if ~exist('lxl','var') || isempty(lxl), 
    lxL = 1:sum(ind);
    lxN = 'Subjects';
else
    lxN = lxl.Name;
    lxL = lxl.Values;
end

%[s_labelh,s_ind] = sort(labelh,'descend');
switch GraphType

    case 0
        % Print title and legend for current subplot
        if h==1
            hx(3) = title(['Out-of-training (OOT) classifications: ' ...
                'means \pm SD'],'FontSize',TitleSize,'FontWeight',TitleWeight);
        end
        pred    = BinClass.mean_predictions(s);
        % Mean predictions with [95%] confidence interval
        errbarCI2 = BinClass.CI2_predictions(s);
        errbarCI1 = BinClass.CI1_predictions(s);
        %         L = zeros(length(pred(ind)),1); U = L;
        %         for i=1:length(pred(ind))
        %             if errbarCI1(ind(i)) < pred(ind(i))
        %                 li = abs(pred(ind(i)) - errbarCI1(ind(i)));
        %                 ui = abs(errbarCI2(ind(i)) - pred(ind(i)));
        %                 if isempty(ui) || isempty(li), continue; end
        %                 U(i) = ui;
        %                 L(i) = li;
        %             else
        %                 ui = abs(pred(ind(i)) - errbarCI1(ind(i)));
        %                 li = abs(errbarCI2(ind(i)) - pred(ind(i)));
        %                 if isempty(ui) || isempty(li), continue; end
        %                 U(i) = ui;
        %                 L(i) = li;
        %             end
        %         end
        % Compute and plot error bars
        L = pred(ind) - errbarCI1(ind);
        U = errbarCI2(ind) - pred(ind);
        errorbar(lxL,pred(ind), L, U,'k', ...
            'LineWidth',ErrorMarkerWidth,'MarkerSize',ErrorMarkerSize, 'LineStyle','none')
    case 1
         if h==1
            hx(3) = title(['Out-of-training (OOT) classifications: ' ...
                'means with 95% CI'],'FontSize',20,'FontWeight','demi');
        end
        pred    = BinClass.mean_predictions(s);
        errbar  = BinClass.std_predictions(s);
        % Mean predictions with standard deviation
        errorbar(lxL,pred(ind),errbar(ind),'k','LineWidth',ErrorMarkerWidth,'MarkerSize',ErrorMarkerSize, 'LineStyle','none')
    case 2
        if h==1
            hx(3) = title(['Out-of-training (OOT) classifications: ' ...
                'group membership probabilities'],'FontSize',TitleSize,'FontWeight',TitleWeight);
        end
        pred    = BinClass.prob_predictions(s,1);
        %b=plot(lxL,pred(ind,1),colpt{groupind(i)},'MarkerSize',16,'LineWidth',3);
    case 4
        if h==1
            hx(3) = title(['Out-of-training (OOT) classifications: ' ...
                'group membership means'],'FontSize',TitleSize,'FontWeight',TitleWeight);
        end
        pred    = BinClass.mean_predictions(s);
end

predh = pred(ind,:);
%% Set GUI data for current binary comparison
xLimitsVec = 1:numel(predh);
xlim([xLimitsVec(1) xLimitsVec(end)]);
pss = cell(1,numel(predh));
for i=1:numel(pss)
    if labelh(i)>0, expgroupi = groupnames{1}; else expgroupi = groupnames{2}; end
    if sign(predh(i))>0, predgroupi = groupnames{1}; else predgroupi = groupnames{2}; end
    pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g'], cases{i}, expgroupi, predgroupi, predh(i));
end
figdata.x = xLimitsVec';
figdata.y = predh;
figdata.patterntext = pss;
guidata(ha, figdata);

%% Mark groups with color
if exist('ind2','var')
    idx1 = false(length(ind1),1); idx1(ind1) = true; idx1 = idx1 & ind;
    idx2 = false(length(ind2),1); idx2(ind2) = true; idx2 = idx2 & ind;
else
    idx1 = false(length(ind1),1); idx1(ind1) = true; idx1 = idx1 & ind;
    idx2 = false(length(~ind1),1); idx2(~ind1) = true; idx2 = idx2 & ind;
end

divx1 = sum(idx1); xl1 = 1:divx1;
divx2 = divx1 + sum(idx2); xl2 = divx1+1:divx2;

if ~isempty(xl1) && ~isempty(xl2)
    switch GraphType
        case {0,1,4}
            if numel(groupind) == 2
                b(1) = plot(lxL(xl1),predh(xl1),colpt{groupind(1)},...
                    'MarkerSize',DataMarkerSize,...
                    'MarkerFaceColor',colptin{groupind(1)},...
                    'LineWidth',DataMarkerWidth,...
                    'LineStyle','none');
                b(2) = plot(lxL(xl2),predh(xl2),colpt{groupind(2)},...
                    'MarkerSize',DataMarkerSize,...
                    'MarkerFaceColor',colptin{groupind(2)},...
                    'LineWidth',DataMarkerWidth,...
                    'LineStyle','none');
            else
                b(1) = plot(lxL(xl1),predh(xl1),colpt{groupind(1)},...
                    'MarkerSize',DataMarkerSize,...
                    'MarkerFaceColor',colptin{groupind(1)},...
                    'LineWidth',DataMarkerWidth,...
                    'LineStyle','none');
                b(2) = plot(lxL(xl2),predh(xl2),'k',...
                    'MarkerSize', DataMarkerSize,...
                    'LineWidth',DataMarkerWidth,...
                    'LineStyle','none');
            end
            xLimits = get(gca,'XLim'); xLimitsVec = xLimits(1):xLimits(2);
            zeroline = zeros(1,numel(xLimitsVec));
            signpred = sign(predh);
        case 2
            predhx = predh(:,1)-0.5;
            if numel(groupind) == 2
                b(1) = plot(lxL(xl1),predh(xl1),colpt{groupind(1)},...
                    'MarkerSize',DataMarkerSize,...
                    'MarkerFaceColor',colptin{groupind(1)},...
                    'LineWidth',DataMarkerWidth,...
                    'LineStyle','none');
                b(2) = plot(lxL(xl2),predh(xl2),colpt{groupind(2)},...
                    'MarkerSize',DataMarkerSize,...
                    'MarkerFaceColor',colptin{groupind(2)},...
                    'LineWidth',DataMarkerWidth,...
                    'LineStyle','none');
            else
                b(1) = plot(lxL(xl1),predh(xl1),colpt{groupind(1)},...
                    'MarkerSize',DataMarkerSize,...
                    'MarkerLineWidth',DataMarkerWidth,...
                    'LineStyle','none');
                b(2) = plot(lxL(xl2),predh(xl2),'k',...
                    'MarkerSize',DataMarkerSize,...
                    'LineWidth',3,...
                    'LineStyle','none');
            end
            xLimits = get(gca,'XLim'); xLimitsVec = xLimits(1):xLimits(2);
            zeroline = ones(1,numel(xLimitsVec))*0.5;
            signpred = sign(predhx);
    end
    plot(xLimitsVec,zeroline,'k--','LineWidth',ZeroLineWidth)

    % Mark misclassifications with respective color

    if size(labelh,1) == 1, labelh = labelh'; end
    err = signpred ~= labelh(ind);
    idx1 = false(length(err),1); 
    idx1(xl1) = true; 
    idx1 = idx1 & err;
    idx2 = false(length(err),1); 
    idx2(xl2) = 1; 
    idx2 = idx2 & err;
    switch GraphType
        case {0,1,4}
            contigmat = BinClass.contigency;
            x=plot(lxL(idx1),predh(idx1),'kx','MarkerSize',DataMissMarkerSize,'LineWidth',DataMissMarkerWidth);
            if isempty(x)
                x=plot(lxL(idx2),predh(idx2),'kx','MarkerSize',DataMissMarkerSize,'LineWidth',DataMissMarkerWidth);
            else
                plot(lxL(idx2),predh(idx2),'kx','MarkerSize',DataMissMarkerSize,'LineWidth',DataMissMarkerWidth);
            end
            % Set descriptions
            switch SVM.prog
                case {'MikSVM','MKLRVM'}
                    algostr = 'RVM probability';
                case 'LIBSVM'
                    algostr = 'SVM score';
                case 'MVTRVR'
                    algostr = 'RVR score';
                case 'LIBLIN'
                    switch SVM.LIBLIN.classifier
                        case {0,6}
                            algostr = 'LR';
                        otherwise
                            algostr = 'SVM';
                    end
                    switch SVM.LIBLIN.b
                        case 0
                            algostr = [algostr ' Score'];
                        case 1
                            algostr = [algostr ' Probability'];
                    end
            end
        case 2
            contigmat = BinClass.prob_contigency;
            x=plot(lxL(idx1),predh(idx1,1),'kx','MarkerSize',DataMissMarkerSize,'LineWidth',DataMissMarkerWidth);
            if isempty(x)
                x = plot(lxL(idx2),predh(idx2,1),'kx','MarkerSize',DataMissMarkerSize,'LineWidth',DataMissMarkerWidth);
            else
                 plot(lxL(idx2),predh(idx2,1),'kx','MarkerSize',DataMissMarkerSize,'LineWidth',DataMissMarkerWidth)
            end
            algostr = 'OOT-Probability';
    end

    if h==nclass,hx(1) = xlabel(lxN); set(hx(1),'FontSize',AxisLabelSize,'FontWeight',AxisLabelWeight),end
    hx(2) = ylabel(algostr); set(hx(2),'FontSize',AxisLabelSize,'FontWeight',AxisLabelWeight)

    legend([b,x],{groupnames{:},'error'},'Location','BestOutside'); legend('boxon')

    % Print classification performance 
    contig = ['Accuracy: '                num2str(contigmat.acc,'%1.1f')];
    contig = char(contig,['Sensitivity: ' num2str(contigmat.sens,'%1.1f')]);
    contig = char(contig,['Specificity: ' num2str(contigmat.spec,'%1.1f')]);
    contig = char(contig,['BAC: '         num2str(contigmat.BAC,'%1.1f')]);
    contig = char(contig,['AUC: '         num2str(contigmat.AUC,'%1.3f')]);
    contig = char(contig,['MCC: '         num2str(contigmat.MCC,'%1.3f')]);
    contig = char(contig,['PPV: '         num2str(contigmat.PPV,'%1.3f')]);
    contig = char(contig,['NPV: '         num2str(contigmat.NPV,'%1.3f')]);
    contig = char(contig,['FPR: '         num2str(contigmat.FPR,'%1.3f')]);
    contig = char(contig,['LR+: '         num2str(contigmat.pLR,'%1.1f')]);
    contig = char(contig,['LR-: '         num2str(contigmat.nLR,'%1.1f')]);
    %xLimits = get(gca,'XLim'); 
    %xlim([xLimits(1) xLimits(2)+1])
    text(-0.15,0.3,contig,'Units','normalized');
end
