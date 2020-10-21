% =========================================================================
% =                          ROC ANALYSIS                                 =
% =========================================================================
function [ hroc, hroc_random] = display_roc(handles, targets, predictions, axeshdl, clafl, linewidth)

if ~exist('axeshdl','var') || isempty(axeshdl), axeshdl = 'axes2'; end 
if ~exist('clafl','var') || isempty(clafl), clafl = true; end 
if ~exist('linewidth','var') || isempty(linewidth), linewidth = 2; end 

%GraphType = get(handles.selYaxis,'Value');
h_class         = get(handles.popupmenu1,'Value');
h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
h_classlist     = get(handles.popupmenu1,'String');

axes(handles.(axeshdl)); 
if clafl, cla; hold on; end

cl = 'k';
if exist('targets', 'var') && ~isempty(targets)
    for i = 1:size(targets,2)
        [X, Y] = perfcurve2(targets(:,i), predictions(:,1), 1);  
        cl = handles.colptin(i,:);
        hroc(i) = plot(handles.(axeshdl),X, Y, 'Color', cl, 'LineWidth', linewidth); 
    end
else
    if strcmp(handles.modeflag,'classification')
        if strcmpi(h_classlist{h_class},'Multi-group classifier')
            for i = 1:handles.ngroups
                if h_onevsall_val == 1
                    cl = handles.colptin(i,:);
                    curclass = i;
                else
                    cl = handles.colptin(h_onevsall_val-1,:);
                    curclass = h_onevsall_val-1;
                end
                X = handles.MultiClass.X{curclass};
                Y = handles.MultiClass.Y{curclass};
                hroc(i) = plot(handles.(axeshdl),X, Y, 'Color', cl, 'LineWidth', linewidth); 
            end
        else
            X = handles.BinClass{h_class}.X;
            Y = handles.BinClass{h_class}.Y;
            hroc = plot(handles.(axeshdl),X, Y, 'Color', cl, 'LineWidth', linewidth); 
        end
    else
        X = handles.Regr.X;
        Y = handles.Regr.Y;
        hroc = plot(handles.(axeshdl),X, Y, 'Color', cl, 'LineWidth', linewidth); 
    end   
        
        %rocout = roc2( [targets(:,i), predictions(:,i),1);
   
end
xlabel(handles.(axeshdl),'False positive rate'); ylabel(handles.(axeshdl),'True positive rate'); 
handles.axes2.XLabel.Color='k'; 
%handles.axes2.XLabel.FontSize=12;
%handles.axes2.YLabel.FontSize=12;
hroc_random = plot(handles.(axeshdl),[0 1],[0 1],'k-.');
