% =========================================================================
% =                        CLASSIFICATION PLOTS                           =
% =========================================================================
function handles = display_classplot(h, handles)

axes(handles.axes1); cla; hold on

switch handles.tglSort.Value
    case 0
        MS = handles.DataMarkerSize; SrtStr = '';
    case 1
        MS = 6; SrtStr = 'Sorted ';
end

%% Display classification plot
% Define X axis data
if ~isfield(handles.BinClass{h},'Xaxis') || isempty(handles.BinClass{h}.Xaxis), 
    lxL = 1:length(handles.BinClass{h}.labelh);
    lxN = [SrtStr 'Subject No.'];
    AltAx = false;
    set(handles.axes1,'Position', handles.axes1pos_orig)
    set(handles.axes38,'Visible','off'); cla(handles.axes38);
else
    lxL = handles.BinClass{h}.Xaxis;
    lxN = handles.XaxisName;
    AltAx = true;
    set(handles.axes1,'Position', handles.axes1pos_alt);
    set(handles.axes38,'Visible','on')
end

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

% Plot errobars if needed
switch GraphType
    case {1,2,3}
       
        predh = handles.BinClass{h}.mean_predictions;
        switch GraphType
            case 2
                % Mean predictions with 95%-CI
                errbarCI2 = handles.BinClass{h}.CI2_predictions;
                errbarCI1 = handles.BinClass{h}.CI1_predictions;
                switch handles.tglSort.Value
                    case 0 
                        L = predh - errbarCI1; U = errbarCI2 - predh;
                        handles.classplot = errorbar(lxL, predh, ...
                                            L, U, ...
                                            'k','LineWidth',handles.ErrorMarkerWidth, ...
                                            'MarkerSize',handles.ErrorMarkerSize, ...
                                            'LineStyle','none');
                    case 1
                        handles.classplot = plotshaded(lxL,[errbarCI1'; errbarCI2'], 'k');
                end
            case 3
                % Mean predictions with standard deviation
                errbar  = handles.BinClass{h}.std_predictions;
                switch handles.tglSort.Value
                    case 0
                        handles.classplot = errorbar(lxL, predh, ...
                                            errbar, ...
                                            'k','LineWidth',handles.ErrorMarkerWidth, ...
                                            'MarkerSize',handles.ErrorMarkerSize, ...
                                            'LineStyle','none');
                    case 1
                         handles.classplot = plotshaded(lxL,[(predh-errbar)'; (predh+errbar)'], 'k');
                end
                       
        end
    case 4
        % Majority voting probabilities
        predh = handles.BinClass{h}.prob_predictions(:,1);
    case 5
        % Cross-CV2 perm majority voting probabilities (95%-CIs)
        predh = handles.BinClass{h}.CV2grid.mean_predictions;
        errbarCI2 = handles.BinClass{h}.CV2grid.CI2_predictions;
        errbarCI1 = handles.BinClass{h}.CV2grid.CI1_predictions;
        switch handles.tglSort.Value
            case 0 
                L = predh - errbarCI1; U = errbarCI2 - predh;
                handles.classplot = errorbar(lxL, predh, ...
                                    L, U, ...
                                    'k','LineWidth',handles.ErrorMarkerWidth, ...
                                    'MarkerSize',handles.ErrorMarkerSize, ...
                                    'LineStyle','none');
            case 1
                handles.classplot = plotshaded(lxL,[errbarCI1'; errbarCI2'], 'k');
        end
    case 6
         % Cross-CV2 perm majority voting probabilities with standard deviation
        predh = handles.BinClass{h}.CV2grid.mean_predictions;
        errbar  = handles.BinClass{h}.CV2grid.std_predictions;
        switch handles.tglSort.Value
            case 0
                handles.classplot = errorbar(lxL, predh, ...
                                    errbar, ...
                                    'k','LineWidth',handles.ErrorMarkerWidth, ...
                                    'MarkerSize',handles.ErrorMarkerSize, ...
                                    'LineStyle','none');
            case 1
                 handles.classplot = plotshaded(lxL,[(predh-errbar)'; (predh+errbar)'], 'k');
        end
end

% Define X axis scaling
if ~AltAx
    r = 1;
    xLimitsVec = 1:numel(predh);
    xLimitsVecInfo = xLimitsVec';
else
    r = (nk_Range(handles.BinClass{h}.Xaxis)/100)*5;
    xLimitsVec = min(handles.BinClass{h}.Xaxis):max(handles.BinClass{h}.Xaxis);
    xLimitsVecInfo = handles.BinClass{h}.Xaxis;
end
XLIMS = [xLimitsVec(1)-r xLimitsVec(end)+r];
xlim(handles.axes1, XLIMS);
ylim(handles.axes1, 'auto');

% Define textbox info data 
pss = cell(1,numel(predh)); psslen=0;
if handles.BinClass{h}.CoxMode
    offs = handles.BinClass{h}.mean_cutoff_probabilities;
else
    offs = zeros(size(predh,1),1);
end

for i=1:numel(pss)
    if handles.BinClass{h}.labelh(i) > 0, 
        expgroupi = handles.BinClass{h}.groupnames{1}; 
    else
        expgroupi = handles.BinClass{h}.groupnames{2};
    end
    if predh(i) > offs(i),
       predgroupi = handles.BinClass{h}.groupnames{1}; 
    else
       predgroupi = handles.BinClass{h}.groupnames{2}; 
    end
    if  handles.BinClass{h}.CoxMode
        pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g' ...
                    '\nCutoff: %g'], handles.BinClass{h}.cases{i}, expgroupi, predgroupi, predh(i), offs(i));
    else
        pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g'], handles.BinClass{h}.cases{i}, expgroupi, predgroupi, predh(i));
    end
    if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.cases       = handles.BinClass{h}.cases;
figdata.x           = xLimitsVecInfo;
figdata.y           = predh;
figdata.patterntext = pss;
figdata.parentui    = handles.pnBinary;
figdata.pnpos       = handles.pnBinary.Position;
figdata.figpos      = handles.figure1.Position;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... %'VerticalAlign', 'Top', ...
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 6, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
                        
set(handles.axes1,'UserData',figdata);

% Mark groups with color
if isfield(handles.BinClass{h},'ind2')
    id1 = handles.BinClass{h}.ind1; 
    id2 = handles.BinClass{h}.ind2; 
else
    id1 = handles.BinClass{h}.ind1; 
    id2 = ~handles.BinClass{h}.ind1; 
end

divx1 = sum(id1); xl1 = 1:divx1;
divx2 = divx1 + sum(id2); xl2 = divx1+1:divx2;

if size(handles.BinClass{h}.labelh,1) == 1, 
    labelh = handles.BinClass{h}.labelh'; 
else
    labelh = handles.BinClass{h}.labelh;
end
if GraphType > 3, 
    signpred = sign(predh-0.5);
else
    signpred = sign(predh-offs);
end
err = signpred ~= labelh;
idx1 = id1 & ~err; idx2 = id2 & ~err; b(1) = 0; b(2)= 0;

if sum(idx1)
    b(1) = plot(lxL(idx1),predh(idx1),handles.colpt,...
        'MarkerSize',MS,...
        'Color',handles.colptin(handles.BinClass{h}.groupind(1),:),...
        'MarkerFaceColor',handles.colptin(handles.BinClass{h}.groupind(1),:),...
        'LineWidth',handles.DataMarkerWidth,...
        'LineStyle','none');    
else
    b(1) = plot(1,NaN,'LineStyle','none');
end

if sum(idx2)
    if ~handles.BinClass{h}.one_vs_all
        CLP = handles.colpt;
        CLR = handles.colptin(handles.BinClass{h}.groupind(2),:);
    else
        CLP = 'o';
        CLR = rgb('DarkGrey');
    end
    b(2) = plot(lxL(idx2),predh(idx2),CLP,...
        'MarkerSize',MS,...
        'Color',CLR,...
        'MarkerFaceColor', CLR,...
        'LineWidth',handles.DataMarkerWidth,...
        'LineStyle','none');
else
    b(2) = plot(1,NaN,'LineStyle','none');
end

if GraphType >3 
    probfx = 0.5;
    predhx = predh-0.5;
    handles.BinClass{h}.predh = predhx;
else
    probfx = 0;
    handles.BinClass{h}.predh = predh;
end
xLimits = get(handles.axes1,'XLim'); xLimitsVec = xLimits(1):xLimits(2);
zeroline = ones(1,numel(xLimitsVec))*probfx;
if handles.BinClass{h}.CoxMode
    zeroline = handles.BinClass{h}.mean_cutoff_probabilities;
    xLimitsVec(1)=[]; xLimitsVec(end)=[];
end
    
ide1 = id1 & err; ide2 = id2 & err;

x1 = plot(lxL(ide1),predh(ide1), '*', 'Color', handles.colptin(handles.BinClass{h}.groupind(1),:),'MarkerSize',handles.DataMissMarkerSize,'LineWidth',handles.DataMissMarkerWidth);
if handles.BinClass{h}.one_vs_all 
    Color2 = rgb('DarkGrey');
else
    Color2 = handles.colptin(handles.BinClass{h}.groupind(2),:);
end
x2 = plot(lxL(ide2),predh(ide2), '*', 'Color', Color2,'MarkerSize',handles.DataMissMarkerSize,'LineWidth',handles.DataMissMarkerWidth);  
handlevec = [b,x1,x2];
legendvec = [handles.BinClass{h}.groupnames(:)',{'misclassified'}, {'misclassified'}];
handles.axes1.XTickMode='auto'; 
handles.axes1.YGrid='off'; 
handles.axes1.XGrid='off'; 

if AltAx,
    % Display regression lines for alternative X Axis
    xl      = get(gca,'Xlim');
    p       = polyfit(lxL(id1),predh(id1),1);     % p returns 2 coefficients fitting r = a_1 * x + a_2
    rho(1)  = corr(lxL(id1),predh(id1));
    r       = p(1) .* xl + p(2);                   % compute a new vector r that has matching datapoints in x
    hl(1)   = plot(xl,r, handles.colpt, 'LineWidth', 3, 'Color', handles.colptin(handles.BinClass{h}.groupind(1),:));
    p       = polyfit(lxL(id2),predh(id2),1);     % p returns 2 coefficients fitting r = a_1 * x + a_2
    rho(2)  = corr(lxL(id2),predh(id2));
    r       = p(1) .* xl + p(2);                   % compute a new vector r that has matching datapoints in x
    hl(2)   = plot(xl,r, handles.colpt, 'LineWidth', 3, 'Color', handles.colptin(handles.BinClass{h}.groupind(2),:));
    handlevec = [handlevec hl];
    legendvec = [legendvec, {['r = ' num2str(rho(1),'%1.2f')]}, {['r = ' num2str(rho(2),'%1.2f')]} ];
    % Display misclassification histogram analysis 
    % Group 1 misclassification histogram
    [err_hist1, Bins1] = hist3([lxL(id1) err(id1)],[10 2]); 
    perr_hist1 = err_hist1(:,2) ./ sum(err_hist1,2);
     % Group 2 misclassification histogram
    [err_hist2, Bins2] = hist3([lxL(id2) err(id2)], Bins1); 
    perr_hist2 = err_hist2(:,2) ./ sum(err_hist2,2);
    axes(handles.axes38); 
    bar(handles.axes38, Bins1{1}, perr_hist1,'BarWidth',1,'FaceColor', handles.colptin(handles.BinClass{h}.groupind(1),:),'FaceAlpha',0.5); 
    xlim(handles.axes38, XLIMS);
    ylim(handles.axes38, [0 1]);
    set(handles.axes38, ...%'FontSize', handles.AxisTickSize, ...
        'FontWeight', handles.AxisTickWeight, ...
        'LineWidth', handles.AxisLineWidth);
    ylabel(handles.axes38,'% Misclassified / Bin');
    hold(handles.axes38,'on'); 
    if handles.BinClass{h}.one_vs_all 
        Color2 = rgb('DarkGrey');
    else
        Color2 = handles.colptin(handles.BinClass{h}.groupind(2),:);
    end
    bar(handles.axes38, Bins1{1}, perr_hist2,'BarWidth',1,'FaceColor', Color2,'FaceAlpha',0.5); 
    hold(handles.axes38,'off'); 
    axes(handles.axes1);
    
    % Kolomogorv-Smirnov-Tests if available
%         if exist('kstest2','file')
%             pnull_hist1 = repmat(sum(err_hist1(:,2))/10,10,1)./sum(err(id1));
%             alt_hist1   = err_hist1(:,2)/sum(err(id1));
%             bs1 = nk_BattacharyyaCoef( alt_hist1, pnull_hist1);
%             pnull_hist2 = repmat(sum(err_hist2(:,2))/10,10,1)./sum(err(id2));
%             alt_hist2   = err_hist2(:,2)/sum(err(id2));
%             bs2 = nk_BattacharyyaCoef( alt_hist2, pnull_hist2);
%             bs3 = nk_BattacharyyaCoef( alt_hist1, alt_hist2);
%             fprintf('\n\n'); cprintf('*black','Battacharyya tests for inequality of distributions');
%             fprintf('\n'); cprintf('*black','========================================================');
%             fprintf('\nGroup %s vs equality:  Battacharyya Coefficient = %1.3f',  handles.BinClass{h}.groupnames{1}, bs1);
%             fprintf('\nGroup %s vs equality:  Battacharyya Coefficient = %1.3f',  handles.BinClass{h}.groupnames{2}, bs2);
%             fprintf('\nGroup %s vs. Group %s: Battacharyya Coefficient = %1.3f', handles.BinClass{h}.groupnames{1}, handles.BinClass{h}.groupnames{2}, bs3);
%             fprintf('\n'); cprintf('*black','========================================================');
%             fprintf('\n')
%         end
end

plot(xLimitsVec,zeroline,'LineWidth',handles.ZeroLineWidth,'Color',rgb('Grey'))
if handles.BinClass{h}.CoxMode
    %plotshaded(lxL,[(zeroline+(offs/2-zeroline))'; (zeroline-(offs/2-zeroline))'], 'k')
end
switch GraphType

    case {1,2,3}
        switch handles.params.TrainParam.SVM.prog
            case {'MikSVM','MKLRVM'}
                algostr = 'RVM probability';
            case 'LIBSVM'
                algostr = 'SVM score';
            case 'MVTRVR'
                algostr = 'RVR score';
            case 'MEXELM';
                algostr = 'ELM score';
            case 'LIBLIN'
                switch handles.params.TrainParam.SVM.LIBLIN.classifier
                    case {0,6}
                        algostr = 'LR';
                    otherwise
                        algostr = 'SVM';
                end
                switch handles.params.TrainParam.SVM.LIBLIN.b
                    case 0
                        algostr = [algostr ' Score'];
                    case 1
                        algostr = [algostr ' Probability'];
                end
            case 'matLRN'
                algostr = sprintf('matLearn [ %s ]',char(handles.params.TrainParam.SVM.matLRN.algo));
            case 'WBLCOX'
                algostr = sprintf('Willbur-Cox proportional hazards model: predicted risk');
            otherwise
                algostr = [handles.params.TrainParam.SVM.prog ' score'];
        end
    case 4
        algostr = 'OOT-Probability';
    case 5
        algostr = 'Mean OOT-Probability (95%-CIs)';
    case 6
        algostr = 'Mean OOT-Probability (SD)';
end
hx(1) = xlabel(lxN); 
set(hx(1), ...%'FontSize',handles.AxisLabelSize-2,
    'FontWeight',handles.AxisLabelWeight);
hx(2) = ylabel([SrtStr algostr]); 
set(hx(2), ...%'FontSize',handles.AxisLabelSize-2, ...
    'FontWeight',handles.AxisLabelWeight);
handles.legend_classplot = legend(handlevec,legendvec, 'Location','Best', 'FontSize', 8,'LineWidth',1);%,'FontSize',handles.LegendFontSize); 
legend('boxon')


%% Display Contigency data
%if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
handles = display_contigmat(handles);

%% Display ROC
if isfield(handles.BinClass{h}.contingency,'AUC')
  flg='on';
    [handles.hroc, handles.hroc_random] = display_roc(handles);
else
  flg = 'off';
end     
handles.axes1.XTickLabelMode = 'auto';
handles.pnRocCmds.Visible = flg;
handles.pnPieCmds.Visible = flg;
handles.pnContigCmds.Visible = flg;
handles.axes20.Visible = flg;
%handles.cmdMetricExport.Visible = flg;
handles.cmdExportAxes20.Visible = flg;

%% Display pie charts
[handles.h1pie, handles.h2pie] = display_piecharts(handles);

%% Display contingency plot
handles.h_contig = display_contigplot(handles, [], handles.BinClass{h}.groupnames);
    
