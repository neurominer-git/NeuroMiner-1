% =========================================================================
% =                            MODELPERF PLOT                             =
% =========================================================================
function handles = display_modelperf(handles)

multfl = false;
predind         = get(handles.popupmenu1,'Value');
h               = predind;
predstr         = get(handles.popupmenu1,'String');
h_list          = get(handles.selModelMeasures,'String');
h_val           = get(handles.selModelMeasures,'Value');
grd             = handles.grid;
meastype        = h_list{h_val};
[~,yl]          = nk_GetScaleYAxisLabel(handles.params.TrainParam.SVM); 
if strcmpi(predstr{predind},'Multi-group classifier')
    h=1; multfl=true;
    switch meastype
        case 'Multi-class cross-validation performance'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam'); 
            P = [ grd.MultiCVPerf grd.MultiTSPerf ]; SEM = [ grd.seMultiCVPerf grd.seMultiTSPerf ];
            if strcmpi(pardesc.GridParam,'BAC'); perfstr = 'Average Binary BAC [%]'; else perfstr = 'Multi-class Accuracy [%]'; end
            ylb = ['Multi-class CV1/CV2-test performance [' perfstr ' ]'];
            lgstr{1} = 'Multi-class CV1 test performance';
            lgstr{2} = 'Multi-class CV2 test performance';
            
        case 'Multi-class generalization error'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
            P = grd.MultiERR_CVTSPerf; SEM = grd.seMultiERR_CVTSPerf; 
            ylb = ['Multi-group generalization error [CV1-test - CV2-test] [' pardesc.GridParam ' ]'];
           
        case 'Multi-class ensemble diversity'
            P = [ grd.MultiCVDiversity grd.MultiTsDiversity]; SEM = [ grd.seMultiCVDiversity grd.seMultiTsDiversity];
            pardesc = 'Entropy-based ensemble diversity';
            ylb = 'Multi-class ensemble entropy';
            
        case 'Multi-class complexity'
            pardesc = 'Average model complexity' ;
            P = grd.MultiComplexity; SEM = grd.seMultiComplexity;
            ylb = 'Average model complexity [%]';
            
        case 'Multi-class model selection frequency'
            pardesc = 'Average parameter selection frequency' ;
            P = grd.MultiSelNodeFreq*100; SEM = [];
            ylb = 'Multi-class parameter selection frequency [%]';
          
    end
else
    switch meastype

        case 'Cross-validation performances'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
            P = [ grd.mean_CVPerf(:,h,:,handles.curlabel) grd.mean_TSPerf(:,h,:,handles.curlabel) ];
            SEM = [ grd.se_CVPerf(:,h,:,handles.curlabel) grd.se_TSPerf(:,h,:,handles.curlabel) ];
            ylb = ['CV1/CV2-test performance [' pardesc.GridParam ' ]'];
            lgstr{1} = 'CV1 test performance';
            lgstr{2} = 'CV2 test performance';
    
        case 'Generalization error'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
            P = grd.mean_Err_CVTSPerf(:,h,:,handles.curlabel);
            SEM = grd.se_Err_CVTSPerf(:,h,:,handles.curlabel);
            ylb = ['Generalization error [CV1-test - CV2-test] [' pardesc.GridParam ' ]'];

        case 'Model complexity'
            pardesc = 'Percentage of observations involved in decision rule [%]';
            P = grd.mean_Complexity(:,h,:,handles.curlabel);
            SEM = grd.se_Complexity(:,h,:,handles.curlabel);
            ylb = 'Model complexity';

        case 'Ensemble diversity'
            P = [ grd.mean_CVDiversity(:,h,:,handles.curlabel) grd.mean_TsDiversity(:,h,:,handles.curlabel)];
            SEM = [ grd.se_CVDiversity(:,h,:,handles.curlabel) grd.se_TsDiversity(:,h,:,handles.curlabel)];
            pardesc = 'Entropy-based ensemble diversity';
            ylb = 'Binary ensemble entropy';
            
        case 'Model selection frequency'
            pardesc = 'Average parameter selection frequency' ;
            P = grd.SelNodeFreq(:,h,:,handles.curlabel);
            SEM = [];
            ylb = 'Parameter selection frequency [%]';
        
        case 'Overall sequence gain'
             pardesc = 'Average prediction performance gain across sequential examination' ;
             P =  grd.mean_SeqGain(:,h,:,handles.curlabel);
             SEM = grd.se_SeqGain(:,h,:,handles.curlabel);
             perf = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
             ylb = ['CV1-test performance gain [ ' perf.GridParam  ' ]'];
             
        case 'Examination frequencies'
            nE = size(handles.params.TrainParam.SVM.SEQOPT.C,2);
            pardesc = 'Percentage of cases of the population at each examination in the sequence';
            P = grd.mean_SeqExamFreq(:,:,h,:,handles.curlabel);
            SEM = grd.se_SeqExamFreq(:,:,h,:,handles.curlabel);
            ylb = ['Examination frequencies [%]'];
            lgstr = cellstr([repmat('Examination frequency: model #', nE ,1) num2str((1:nE)')]);
            
        case 'Case propagation thresholds'
            nE = size(handles.params.TrainParam.SVM.SEQOPT.C,2)-1;
            pardesc = 'Percentage of cases of the population at each examination in the sequence';
            P = [grd.mean_SeqPercUpper(:,:,h,:,handles.curlabel) -1*grd.mean_SeqPercLower(:,:,h,:,handles.curlabel)];
            SEM = [grd.se_SeqPercUpper(:,:,h,:,handles.curlabel) -1*grd.se_SeqPercLower(:,:,h,:,handles.curlabel)];
            ylb = ['Examination frequencies [%]'];
            lgstr = [repmat('+1 * Upper propagation threshold: model #', nE ,1) num2str((1:nE)')];
            lgstr = [lgstr; repmat('-1 * Lower propagation threshold: model #', nE ,1) num2str((1:nE)')];
            lgstr = cellstr(lgstr);
    end
    
end
handles.currmeas = P;
handles.currmeasdesc = pardesc;

[handles, contfl] = display_SubParam(handles ,'display_modelperf');
if ~contfl, 
    return; 
end

Ps              = handles.MLparams.ParamCombs;
Ps_desc         = handles.MLparams.ParamDesc;
%NumParamDims    = handles.MLparams.NumParamDims;

if isfield(handles,'BinClass') || isfield(handles,'MultiClass')
    if multfl 
        groupdesc = 'Multi-group classification model';
    else
        groupdesc = ['Classification model (' handles.BinClass{h}.description '):'];
    end
elseif isfield(handles,'Regr')
    groupdesc = 'Regression model';
end

[lPi,nPi] = size(P);
if nPi > 1,
    cpt = handles.colptin;
else
    cpt = rgb('Lawngreen');
end

if lPi > 500
    MrkSize = 2; LinSize = 0.1;
elseif lPi > 250
    MrkSize = 5; LinSize = 0.25;
elseif lPi > 150
    MrkSize = 7.5; LinSize = 0.5;
else
    MrkSize = 10; LinSize = 1;
end
   
axes(handles.axes17); cla; hold on
legend('hide')
mlx = Inf; mux = 0;

% Plot performance graphs
switch meastype
    
    case 'Case propagation thresholds'
        x = repmat((0.5:lPi+0.5)',nPi,1); y = P(:); hi = zeros(nPi,1);
        xi = (1:lPi)'; mlx = -100; mux = 100;
        for i=1:2:nPi
             yi = P(:,i:i+1); 
             hi(i:i+1) = bar(xi,yi,'grouped','LineWidth', 1.0);
        end
        % Define x and y axis label
        strp = [];
        for i = 1 : numel(Ps_desc{h})
            strp = sprintf('%s, %s', strp,Ps_desc{h}{i} ); 
        end
        strp = strp(3:end);
        xlabel(['ML parameters combinations: [' strp ']']);
        ylabel(ylb);

    otherwise
        
        x = repmat((1:lPi)',nPi,1); y = P(:); hi = zeros(nPi,1);
        for i=1:nPi
            xi = (1:lPi)'; y2 = P(:,i);
            if exist('SEM','var') && ~isempty(SEM)
                y1 = P(:,i)-SEM(:,i)/2; y3 = (P(:,i)+SEM(:,i)/2);  y1(isnan(y1))=0; y3(isnan(y3))=0;
                mlxi = min(y1); muxi = max(y3);
                hi(i) = plotshaded(xi',[y1'; y2'; y3'], cpt(i,:), MrkSize, LinSize) ;
            else
                hi(i) = plot(xi,y2,'ko-','MarkerFaceColor', cpt(i,:),'MarkerSize', MrkSize, 'MarkerEdgeColor', rgb('Black'), 'LineWidth', 1.0);
                mlxi = nm_nanmin(y2); muxi = nm_nanmax(y2);
            end
            if mlxi < mlx, mlx = mlxi; end
            if muxi > mux, mux = muxi; end    
        end
        
        % Set GUI data
        pss = cell(1,numel(y));
        ll = 1;
        for k=1:size(P,2)
            for i=1:size(P,1)
                psj=sprintf(' %s: %1.3f', ylb, P(i,k));
                for j = 1:size(Ps_desc{h},2)
                    if iscell(Ps{h}(i,j)),
                        ijPs = Ps{h}{i,j};
                    else
                        ijPs = Ps{h}(i,j);
                    end
                    if ~ischar(ijPs), ijPs = num2str(ijPs); end
                    psj = sprintf('%s\n%s: %s',psj, Ps_desc{h}{j}, ijPs);
                end
                psj = sprintf('Parameter combination %g\n%s',ll, psj(2:end));
                pss{ll} = sprintf('%s',psj);
                ll=ll+1;
            end
        end

        hText = uicontrol('Style','text', 'String', pss{1}, 'FontSize',11, 'Units','normalized', 'Parent',gcf, 'Visible','off'); 
        figdata.x = x;
        figdata.y = y;
        figdata.patterntext = repmat(pss', nPi,1);
        figdata.parentui    = handles.pnModelPerf;
        figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent+0.01, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
        figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                                    'Interpreter','none', ... 
                                    'Color', 'black', ...
                                    'BackgroundColor',[.6 .7 .6], ...
                                    'Position', [0 0 0.99 0.99], ...
                                    'EdgeColor', [.6 .7 .6], ...
                                    'LineWidth', 0.1, ...
                                    'Margin', 5, ...
                                    'FitBoxToText','on', ...
                                    'Visible','off');

        set(handles.axes17,'UserData',figdata);

        % Define x and y axis label
        strp = [];
        for i = 1 : numel(Ps_desc{h})
            strp = sprintf('%s, %s', strp,Ps_desc{h}{i} ); 
        end
        strp = strp(3:end);
        xlabel(['ML parameters combinations: [' strp ']']);
        ylabel(ylb);

        % Mark maximum and minimum values
        for i = 1 : nPi
            yi = P(:,i);
            [maxY, maxX] = max(yi);
            [minY, minX] = min(yi);
            plot(maxX,maxY,'ro','MarkerSize',MrkSize+5, 'LineWidth', 3, 'Parent', handles.axes17)
            plot(minX,minY,'bo','MarkerSize',MrkSize+5, 'LineWidth', 3, 'Parent', handles.axes17)
        end
        
end

% Plot markers
%plot(x, y, 'k.', 'MarkerSize', 1

% Scale axes
if numel(unique(x)) > 1
    xlim([min(x) max(x)])
end
if numel(unique(y)) > 1,
    ylim([mlx mux])
end
handles.axes17.XTickMode = 'auto'; handles.axes17.XTickLabelMode = 'auto';

switch handles.NM.modeflag
    case 'classification'
        if multfl
            ct = handles.MultiClass;
        else
            ct = handles.BinClass{h};
        end
    case 'regression'
        ct = handles.Regr;
end

axes(handles.axes35); cla; colorbar('delete')
axes(handles.axes37); cla; colorbar('delete')

S = handles.axes17.Position;
T_lower = handles.axes37.Position;
T_upper = handles.axes35.Position;

h_sub = S(4)/2-0.03;

handles.axes37.Position = [T_lower(1) T_lower(2) T_lower(3) h_sub];
handles.axes35.Position = [T_upper(1) T_lower(2)+h_sub+0.06 T_upper(3) h_sub];
handles.axes35.XLabel.Visible='off';

switch meastype

    case {'Cross-validation performances', 'Multi-class cross-validation performance'}
        handles.axes17.Position(3) = handles.axes37.Position(1) - handles.axes17.Position(1) - 0.05 ;
        handles.axes35.Visible = 'on'; handles.axes35.YLabel.Visible = 'on'; handles.axes35.XLabel.Visible = 'off'; handles.axes35.Title.Visible = 'on'; 
        handles.axes37.Visible = 'on'; handles.axes37.YLabel.Visible = 'on'; handles.axes37.XLabel.Visible = 'on'; handles.axes37.Title.Visible = 'on'; 
        handles.cmdExportModelPerfC1.Visible='on'; handles.cmdExportModelPerfC2.Visible='on';
        % Print CV2 partition data
        
        axes(handles.axes35); colormap(jet); hold on
        mux = max(ct.best_TS(:)); mlx = min(ct.best_TS(:));
        imagesc(ct.best_TR); 
        ylim([0.5 size(ct.best_TR,1)+0.5]); xlim([0.5 size(ct.best_TR,2)+0.5]); caxis([mlx mux]);
        handles.axes35.YTick = 0:size(ct.best_TR,1);
        title('Mean CV1 Performance at CV2 partition level')
        ylabel('Permutations'); xlabel('Folds')
        
        axes(handles.axes37); colormap(jet); hold on
        imagesc(ct.best_TS);
        ylim([0.5 size(ct.best_TR,1)+0.5]); xlim([0.5 size(ct.best_TR,2)+0.5]); caxis([mlx mux]);
        handles.axes37.YTick = 0:size(ct.best_TR,1);
        hc = colorbar(handles.axes37,'location','SouthOutside');
        ylabel('Permutations'); xlabel('Folds')
        title('Mean CV2 Performance at CV2 partition level')
        hc.Position(2) = 0.04; hc.Position(4) = .015;
        
    case {'Model complexity','Generalization error', 'Multi-class model complexity', 'Multi-class generalization error'}
        handles.axes17.Position(3) = handles.axes37.Position(1) - handles.axes17.Position(1) - 0.05 ;
        handles.axes35.Visible = 'off'; handles.axes35.YLabel.Visible = 'off'; handles.axes35.XLabel.Visible = 'off'; handles.axes35.Title.Visible = 'off'; 
        handles.axes37.Visible = 'on'; handles.axes37.YLabel.Visible = 'on'; handles.axes37.XLabel.Visible = 'on'; handles.axes37.Title.Visible = 'on'; 
        handles.cmdExportModelPerfC1.Visible='off'; handles.cmdExportModelPerfC2.Visible='on';
        
        axes(handles.axes37); 
        colormap(jet); hold on
        switch meastype
            case {'Model complexity','Multi-class model complexity'}
                imagesc(ct.best_Complexity);
                caxis([min(ct.best_Complexity(:)) max(ct.best_Complexity(:))]);
                if handles.params.TrainParam.GRD.OptRegul.flag 
                    title('Mean regularization function at CV2 partition level')
                else
                    title('Mean model complexity at CV2 partition level')
                end
            case {'Generalization error', 'Multi-class generalization error'}
                imagesc(ct.best_Error);
                title('Mean generalization error at CV2 partition level')
                caxis([min(ct.best_Error(:)) max(ct.best_Error(:))]);
        end
        ylim([0.5 size(ct.best_Complexity,1)+0.5]); xlim([0.5 size(ct.best_Complexity,2)+0.5]);
        handles.axes37.YTick = 0:size(ct.best_TR,1);
        ylabel('Permutations'); xlabel('Folds')
        hc = colorbar(handles.axes37,'location','SouthOutside');
        sz1 = handles.axes35.Position; sz2 = handles.axes37.Position;
        handles.axes37.Position = [sz1(1) sz2(2) sz1(3) sz2(4)];
        hc.Position(2) = 0.04; hc.Position(4) = .015;
        
    otherwise
        handles.axes17.Position(3) = handles.axes37.Position(1)+ handles.axes37.Position(3) - handles.axes17.Position(1) ;
        handles.axes35.Visible = 'off'; handles.axes35.YLabel.Visible = 'off'; handles.axes35.XLabel.Visible = 'off'; handles.axes35.Title.Visible = 'off'; 
        handles.axes37.Visible = 'off'; handles.axes37.YLabel.Visible = 'off'; handles.axes37.XLabel.Visible = 'off'; handles.axes37.Title.Visible = 'off'; 
        handles.cmdExportModelPerfC1.Visible='off'; handles.cmdExportModelPerfC2.Visible='off';
        
end
handles.cmdExportModelPerfVec.Position(1) = handles.axes17.Position(1)+handles.axes17.Position(3)+0.005;
if exist('lgstr','var') && ~isempty(lgstr), 
    handles.legend_modelperf = legend(hi,lgstr,'Location','Best','FontSize',handles.LegendFontSize,'LineWidth',1); 
else
    legend('hide')
    handles.legend_modelperf = [];
end