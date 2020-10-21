% =========================================================================
% =                             REGRESSION PLOT                           =
% =========================================================================
function handles = display_regrplot(handles, MarkFlag)

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');
label     = handles.Regr.labels;
if ~isfield(handles.Regr,'Xaxis') || isempty(handles.Regr.Xaxis),
    indnan    = ~isnan(label);
    label     = label(indnan);
    lxL = label;
    lxN = 'Observed targets';
else
    indnan  = ~isnan(handles.Regr.Xaxis);
    label  = label(indnan);
    lxL = handles.Regr.Xaxis(indnan);
    lxN = handles.XaxisName;
end

if ~exist('MarkFlag','var') || isempty(MarkFlag), MarkFlag=false; end

axes(handles.axes1);
uistack(handles.axes1,'top')
cla;
hold on

pred    = handles.Regr.mean_predictions(indnan);
errbar  = handles.Regr.std_predictions(indnan);
ind     = handles.Regr.index_predictions(indnan);
lx      = length(label);

if isfield(handles.Regr,'grouping') && MarkFlag
    ngroups = length(unique(handles.Regr.grouping));
    grouping = handles.Regr.grouping;
    markgroups = true;
else
    ngroups = 1; grouping = ones(lx,1); markgroups = false;
end

% Define textbox info data 
findnan = find(indnan);
pss = cell(1,numel(findnan)); psslen=0;
for i=1:numel(findnan)
      pss{i} = sprintf(['Subject ID: %s' ...
            '\nObserved Target: %g' ...
            '\nPredicted Target: %g\n'], handles.subjects{findnan(i)}, label(i), pred(i));
     if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.x = lxL;
figdata.y = pred;
figdata.patterntext = pss;
figdata.parentui = handles.pnBinary;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... %'VerticalAlign', 'Top', ...
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);

switch GraphType
    
    case 1
        lgstr{1} = '$\mathbf{Median_{OOT}}$';

    case 2
        lgstr{1} = '$\mathbf{Median_{OOT} \pm 95\%CI}$';
        % Mean predictions with [95%] confidence interval
        errbarCI2 = handles.Regr.CI2_predictions(indnan);
        errbarCI1 = handles.Regr.CI1_predictions(indnan);
        L = pred - errbarCI1; U = errbarCI2 - pred;
        handles.he  = errorbar(lxL(ind),pred(ind), L(ind), U(ind),'ko','LineWidth',0.5,'MarkerSize',9);
    
    case 3
        lgstr{1} = '$\mathbf{Median_{OOT} \pm SD}$';
        % Mean predictions with standard deviation
        handles.he = errorbar(lxL(ind),pred(ind),errbar(ind),'ko','LineWidth',0.5,'MarkerSize',9);
end

handles.regrplot = plot(lxL(ind),pred(ind),'o','MarkerFaceColor','b','MarkerSize',9);

% Mark points according to "grouping"
handles.hg = [];
if markgroups
    for i = 1:ngroups
        indg = (ind == 1 & grouping == i);
        handles.hg(i) = plot(lxL(indg),pred(indg),handles.colpt{i},'LineWidth',1,'MarkerSize',9);
    end
    if isfield(handles.Regr,'grouping') && ngroups>1 
        for i=1:ngroups
            lgstr{i+1} = handles.Regr.groupnames{i};
        end
    end
end

xstep = nk_Range(lxL)/100;  r = (xstep)*5; 
xLimitsVec = min(lxL):xstep:max(lxL);
xlim([xLimitsVec(1)-r xLimitsVec(end)+r]);
ylim('auto');

% Add regression line to plot
xy = min(lxL(ind)):(max(lxL(ind))-min(lxL(ind)))/sum(ind):max(lxL(ind));
try % Statistics toolbox available
    [p,s] = polyfit(lxL(ind),pred(ind),1);
    [yhat,dy] = polyconf(p,xy,s,'predopt','curve');
    handles.hline = plot(xy,yhat,'k-','LineWidth',2);
    handles.hline_CI = plotshaded(xy,[yhat+dy; yhat-dy],'b');
    lgstr{end+1} = '$\mathbf{\hat{y}_{linear}}$';
    lgstr{end+1} = '$\mathbf{\hat{y}_{95\%CI}}$';
catch % or not
    handles.hline = lsline;
    handles.hline_CI = [];
    lgstr{end+1} = '$\mathbf{\hat{y}_{linear}}$';
end
% Prepare legend
switch GraphType
    case 1
        hdlvec = [handles.regrplot handles.hline handles.hline_CI ];
    case {2,3}
        hdlvec = [handles.he handles.hg handles.hline handles.hline_CI];
end
handles.legend_classplot = legend(hdlvec, lgstr, 'Location','southeast', 'LineWidth', 1, 'Interpreter','latex');

% Set descriptions
xlabel(lxN)
ylabel('Predicted targets')

axes(handles.axes5)
if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
    
% Binarize at median
m = nm_nanmedian(label); set(handles.txtBinarize,'String',m);

