function handles = display_tablegenerator(handles)
 
if ~isfield(handles,'PerfTabWin')
    WinTag = 'PerfTabWin';
    h = findobj('Tag',WinTag);
    if isempty(h)
       handles = init_fig(handles, WinTag);
    else
       handles.PerfTab.Win = h; 
    end
end
h = handles.PerfTab.Win ;
set(0,'CurrentFigure',h)

function handles = init_fig(handles, WinTag)

handles.PerfTab.Win = figure('NumberTitle','off','Name','NM Performance Tabulator', 'Tag' ,WinTag,'MenuBar','none');

cnt=0; analyses=[];

for i=1:numel(handles.NM.analysis)
    if handles.NM.analysis{i}.status
        cnt = cnt+1;
        %analyses{cnt} = handles.NM.analysis{i}.id;
        analysescnt{cnt} = num2str(cnt);
        data{cnt, 1} = handles.NM.analysis{i}.id;
        data{cnt, 2} = '';
        data{cnt, 3} = false;
    end
end

switch handles.modeflag
    case 'classification'
        mfld = 'BinClass';
    case 'regression'
        mfld = 'Regr';
end

fld = fieldnames(handles.NM.analysis{1}.GDdims{1}.(mfld){1}.contigency);
fcnt = 0;
for i=1:numel(fld)
    if strcmp(fld{i},'X'), 
        continue, 
    end
    fcnt=fcnt+1;
    fldcnt{fcnt} = num2str(fcnt);
    flddata{fcnt,1} = fld{i}; flddata{fcnt,2} = false; 
end
handles.PerfTab.analysisselect = uitable('units', 'normalized', 'position', [0.05, 0.2, 0.9, 0.425], ...
                                        'ColumnName', {'Analyses','Alias','Select'}, ... 
                                        'ColumnFormat',{'char','char','logical'},...
                                        'ColumnEditable', [false, true, true],...
                                        'ColumnWidth',{300, 100, 'auto'}, ...
                                        'RowName', analyses,...
                                        'data', data);
handles.PerfTab.perfselect = uitable('units', 'normalized', 'position', [0.05, 0.65, 0.9, 0.30], ...
                                        'ColumnName', {'Metric','Select'}, ... 
                                        'ColumnFormat',{'char','logical'},...
                                        'ColumnEditable', [false, true],...
                                        'RowName', fldcnt,...
                                        'data', flddata);
handles.PerfTab.subgroupvartext = uicontrol('Style','edit', ...
                                        'units','normalized', ...
                                        'Position',[0.05 0.12 0.278 0.06], 'ToolTip','enter name of categorical variable for subgroup analysis [numeric array or cell array of strings');
handles.PerfTab.subgroupvarnamestext = uicontrol('Style','edit', ...
                                        'units','normalized', ...
                                        'Position',[0.34 0.12 0.288 0.06], 'ToolTip','enter name of cell array variable or provide comma delimited strings describing categories of numeric subgroup array');
handles.PerfTab.subgroupvarclear     = uicontrol('Style','pushbutton', ...
                                        'units','normalized', ...
                                        'Position',[0.64 0.12 0.3122 0.067], ...
                                        'String','Clear subgroup settings', ...
                                        'FontWeight', 'normal',...
                                        'BackgroundColor', rgb('OldLace'), ...
                                        'Callback', {@clearsubgroupsettings, handles});   
handles.PerfTab.fileseltext = uicontrol('Style','edit', ...
                                        'units','normalized', ...
                                        'Position',[0.05 0.04 0.58 0.06], 'ToolTip', 'enter path and name of table file to be generated');
handles.PerfTab.fileseldlg  = uicontrol('Style','pushbutton', ...
                                        'units','normalized',...
                                        'Position',[0.64 0.04 0.15 0.067], ...
                                        'String','Save as', ...
                                        'Callback', {@saveas,handles}, 'ToolTip', 'open file dialog to generate path and name of table file to be generated');
handles.PerfTab.tabulate    = uicontrol('Style','pushbutton', ...
                                        'units','normalized', ...
                                        'Position',[0.80 0.04 0.15 0.067], ...
                                        'String','Tabulate', ...
                                        'FontWeight', 'bold',...
                                        'BackgroundColor', rgb('lightblue'), ...
                                        'Callback', {@tabulate, handles});

function handles = saveas(src, evt, handles)

if ispc 
    ext = '*.xlsx';
else
    ext = '*.csv';
end

[FileName,PathName] = uiputfile(ext,'Save performance table','PerfTable');
handles.PerfTab.fileseltext.String = fullfile(PathName, FileName);

function handles = clearsubgroupsettings(src, evt, handles)


function handles = tabulate(src, evt, handles)

% Check whether path can be created
pth = fileparts(handles.PerfTab.fileseltext.String);
if isempty(handles.PerfTab.fileseltext.String) || isempty(pth), errordlg('Provide a valid output path before tabulating the data.'); return; end

if ~isfolder(pth)
    [status, msg] = mkdir(pth);
    if status
        rmdir(pth)
    else
        errordlg(msg)
        return
    end
end

switch handles.modeflag
    case 'classification'
        mfld = 'BinClass';
    case 'regression'
        mfld = 'Regr';
end
% Retrieve user input
AnalysisSelection = cell2mat(handles.PerfTab.analysisselect.Data(:,3));
MetricSelection = handles.PerfTab.perfselect.Data(cell2mat(handles.PerfTab.perfselect.Data(:,2)),1);
AnalysisStrings = handles.PerfTab.analysisselect.Data(:,1);
AnalysisAliasStrings = handles.PerfTab.analysisselect.Data(:,2);
I = strcmp(AnalysisAliasStrings,'');
AnalysisAliasStrings(I) = AnalysisStrings(I);
AnalysisAliasStringsSel = AnalysisAliasStrings(AnalysisSelection);
SubgroupSelectionVarName = handles.PerfTab.subgroupvartext.String;
SubgroupSelection=[]; strmode = false;
% Check user input!
if ~any(AnalysisSelection)
    errordlg('You have to select at least one analysis from the list')
    return
elseif isempty(MetricSelection), 
    errordlg('You have to select at least one performance metric from the list')
    return; 
elseif ~isempty(SubgroupSelectionVarName) || ~strcmp(SubgroupSelectionVarName,'')
    %% Subgroup analysis preparations
    try
        SubgroupSelection = evalin('base',SubgroupSelectionVarName);
        uSG = unique(SubgroupSelection); 
        if iscell(uSG)
            fNan = strcmp(uSG,''); uSG(fNan)=[];
        else
            uSG(isnan(uSG))=[]; 
        end
        nuSG = numel(uSG);
        if ~isempty(SubgroupSelection)
            if size(SubgroupSelection,1) ~= size(handles.NM.label,1),
                errordlg(sprintf('%s has not the same number of values [%g] as cases in the NM structure!',  SubgroupSelectionVarName, size(handles.NM.label,1)));
                return
            elseif ~isnumeric(SubgroupSelection) && ~iscellstr(SubgroupSelection)
                errordlg(sprintf('%s has to be either as numeric array or a cell array of strings/string array!',  SubgroupSelectionVarName)) 
                return
            elseif isnumeric(SubgroupSelection) 
                SubgroupSelectionNames = handles.PerfTab.subgroupvarnamestext.String;
                if ~isempty(SubgroupSelectionNames ) || ~strcmp(SubgroupSelectionNames ,'')
                    if any(strfind(SubgroupSelectionNames,','))
                        SubgroupSelectionDesc = strsplit(SubgroupSelectionNames,',');
                    elseif any(strfind(SubgroupSelectionNames,'"'))
                        SubgroupSelectionDesc = evalin('base',SubgroupSelectionNames);
                    else
                        SubgroupSelectionDesc = {SubgroupSelectionNames};
                    end
                    nuSGD = numel(SubgroupSelectionDesc);
                    if nuSGD ~= nuSG
                        errordlg(sprintf('The subgroup description variable [%g] must have the same number of entries as the subgroup vector has unique indices [%g]!', nuSGD, nuSG));
                        return
                    end
                end
            elseif iscellstr(SubgroupSelection)
                SubgroupSelectionDesc = uSG;
                strmode = true;
            end
        else
            errordlg(sprintf('%s is empty!',  SubgroupSelectionVarName));
            return
        end
    catch
        errordlg(sprintf('Could not read %s from the MATLAB workspace',  SubgroupSelectionVarName));
        return
    end
end
data_table = cell(1,numel(MetricSelection)+1);
data_table(1,:) = ['Predictors' MetricSelection' ];
cnt=1; fcnt=0; acnt=0;
for i=1:numel(handles.NM.analysis)
    if handles.NM.analysis{i}.status, fcnt=fcnt+1; else; continue; end
    multiflag = false; if isfield(handles.NM.analysis{i}.GDdims{1},'MultiClass'), multiflag=true; end
    if AnalysisSelection(fcnt)
        acnt=acnt+1;
        if i > 1 && multiflag, cnt=cnt+2; else, cnt=cnt+1; end
        data_table{cnt,1} = AnalysisAliasStringsSel{acnt};
        nGD = numel(handles.NM.analysis{i}.GDdims);
        for j = 1:nGD
            if nGD>1, 
                cnt=cnt+1; data_table{cnt,1}=sprintf('Modality #%g', handles.NM.analysis{i}.params.TrainParam.FUSION.M(j)); 
            end
            sfld = handles.NM.analysis{i}.GDdims{j}.(mfld); 
            if multiflag, cnt=cnt+1; data_table{cnt,1} = 'Binary classification results'; end
            for k=1:numel(sfld)
                if multiflag, 
                    cnt=cnt+1; data_table{cnt,1} = sprintf('%s', handles.NM.analysis{i}.params.cv.class{1,1}{k}.groupdesc); 
                end
                if ~isempty(SubgroupSelection)
                    if isfield(sfld{k},'globalcutoff_probabilities')
                        P = sfld{k}.prob_predictions(:,1) - sfld{k}.globalcutoff_probabilities;
                    else
                        P = sfld{k}.mean_predictions;
                    end
                    S = SubgroupSelection(sfld{k}.index_predictions);
                    L = handles.NM.label(sfld{k}.index_predictions);
                    Lx = zeros(size(L));
                    uL = unique(L); vec = [1 -1];
                    for q=1:numel(uL)
                        indq = L == uL(q);
                        Lx(indq) = vec(q);
                    end
                    hdr = data_table{cnt,1};
                    for p=1:nuSG
                        data_table{cnt,1} = sprintf('%s [group=%s]', hdr, SubgroupSelectionDesc{p});
                        if ~strmode
                            indp = S == uSG(p);
                        else
                            indp = strcmp(S,uSG{p});
                        end
                        ALL = ALLPARAM(Lx(indp), P(indp));
                        for l=1:numel(MetricSelection), 
                            data_table{cnt,l+1} = ALL.(MetricSelection{l});
                        end
                        cnt=cnt+1;
                    end
                    cnt=cnt-1;
                else
                     for l=1:numel(MetricSelection), 
                        data_table{cnt,l+1} = sfld{k}.contigency.(MetricSelection{l});
                     end
                end
            end
            if multiflag
                mcfld = handles.NM.analysis{i}.GDdims{j}.MultiClass; 
                cnt=cnt+1; data_table{cnt,1} = 'Multi-group classification results';
                for k=1:numel(sfld)
                    cnt=cnt+1; data_table{cnt,1} = sprintf('Multi-group: %s vs. REST', handles.NM.groupnames{k});
                    for l=1:numel(MetricSelection), data_table{cnt,l+1} = mcfld.class{k}.(MetricSelection{l}); end
                end
            end
        end
    end
end
tbl.colnames = data_table(1,1:end);
tbl.rownames = data_table(2:end,1);
tbl.array    = data_table(2:end,2:end);
[ERR, STATUS, fil, typ] = tbl2file(tbl, handles.PerfTab.fileseltext.String, 'Performance Table');

if STATUS
    msgbox(['Data successfully exported to file ' fil]);
else
    warndlg(sprintf('Data NOT successfully exported to file.\nError: %s', ERR.id));
end
