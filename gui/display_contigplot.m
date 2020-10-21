% =========================================================================
% =                      CONTINGENCY MATRIX PLOT                          =
% =========================================================================
function h = display_contigplot(handles, confmatrix, groupnames)

h_class         = get(handles.popupmenu1,'Value');
h_classlist     = get(handles.popupmenu1,'String');
h_list          = get(handles.selModelMeasures,'String');
h_val           = get(handles.selModelMeasures,'Value');
h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
if ~exist('groupnames','var') || isempty(groupnames), groupnames = handles.NM.groupnames; end
switch h_list{h_val}
    case 'Classification plot'
        if strcmpi(h_classlist{h_class},'Multi-group classifier')
            if h_onevsall_val >1
                if ~exist('confmatrix','var') || isempty(confmatrix)
                    c = handles.MultiClass.class{h_onevsall_val-1};
                    confmatrix = [[c.TP c.FN]; [c.FP c.TN]];
                end
                groupnames = {handles.NM.groupnames{h_onevsall_val-1} , 'REST' };
            else
                if ~exist('confmatrix','var') || isempty(confmatrix)
                    confmatrix = handles.MultiClass.confusion_matrix;
                end
            end
        else
            if ~exist('confmatrix','var') || isempty(confmatrix)
                c = handles.BinClass{h_class}.contingency;
                confmatrix = [[c.TP c.FN]; [c.FP c.TN]];
            end
        end
    case 'Regression plot'
       if ~exist('confmatrix','var') || isempty(confmatrix), confmatrix = handles.Regr.contigmat; end
       thresh = get(handles.txtBinarize,'String');
       groupnames{1} = sprintf('Group 1 (>=%s)',thresh);
       groupnames{2} = sprintf('Group 2 (< %s)',thresh);  
end
n_mat = bsxfun(@rdivide, confmatrix, sum(confmatrix,2))*100;
axes(handles.axes20); h = nk_PlotConfusionMat(n_mat, groupnames, gca, handles.txtPerf);
handles.cmdExportAxes20.Position(2) = handles.axes20.Position(2)+ handles.axes20.Position(4)- handles.cmdExportAxes20.Position(4); 
