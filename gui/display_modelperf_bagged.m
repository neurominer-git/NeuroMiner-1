% =========================================================================
% =                   MODELPERF PLOT (BAGGED ENSEMBLES)                   =
% =========================================================================
function handles = display_modelperf_bagged(handles)

multfl = false;
predind         = get(handles.popupmenu1,'Value');
h               = predind;
predstr         = get(handles.popupmenu1,'String');
h_list          = get(handles.selModelMeasures,'String');
h_val           = get(handles.selModelMeasures,'Value');
meastype        = h_list{h_val};
grd             = handles.GDdims;
axes(handles.axes37); colorbar('delete');

if strcmpi(predstr{predind},'Multi-group classifier')
    h=1; multfl=true;
    switch meastype
        case 'Multi-class cross-validation performance summary'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam'); 
            P = [ grd.best_MultiCVperf_unimodal grd.best_MultiTSperf_unimodal ] ; 
            %SEM = [ grd.best_sdMultiCVperf_unimodal grd.best_sdMultiTSperf_unimodal ];
            if strcmpi(pardesc.GridParam,'BAC'); perfstr = 'Average Binary BAC [%]'; else perfstr = 'Multi-class Accuracy [%]'; end
            ylb = ['Multi-class CV1/CV2 performance [' perfstr ' ]'];
            lgstr{1} = 'Multi-class CV1 test performance';
            lgstr{2} = 'Multi-class CV2 validation performance';
            lgstr{3} = 'Multi-class CV2 validation performance (bagged predictor)';
            
        case 'Multi-class generalization error summary'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
            P = grd.best_sdMultiCVperf_unimodal - grd.best_sdMultiTSperf_unimodal ; 
            %SEM = grd.best_sdMultiCVperf_unimodal - grd.best_sdMultiTSperf_unimodal ;
            ylb = ['Multi-group generalization error [ CV1 - CV2 ] [' pardesc.GridParam ' ]'];
          
    end
else
    switch meastype

        case 'Cross-validation performance summary'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
            P = [ grd.best_CVperf_unimodal(:,h,handles.curlabel) grd.best_TSperf_unimodal(:,h,handles.curlabel) ];
            %SEM = [ grd.best_sdCVperf_unimodal(:,h,handles.curlabel) grd.best_sdCVperf_unimodal(:,h,handles.curlabel) ];
            ylb = ['CV1/CV2 performance [' pardesc.GridParam ' ]'];
            lgstr{1} = 'CV1 test performance (unimodal)';
            lgstr{2} = 'CV2 validation performance (unimodal)';
            lgstr{3} = 'CV2 validation performance (bagged predictor)';
    
        case 'Generalization error summary'
            pardesc = nk_GetParamDescription2([],handles.params.TrainParam,'GridParam');
            P = grd.best_Error_unimodal(:,h,handles.curlabel) ;
            %SEM = grd.best_sdCVperf_unimodal(:,h,handles.curlabel) - grd.best_sdCVperf_unimodal(:,h,handles.curlabel) ;
            ylb = ['Generalization error [ CV1 - CV2 ] [' pardesc.GridParam ' ]'];

        case 'Model complexity summary'
            pardesc = 'Percentage of observations involved in decision rule [%]';
            P = grd.best_Complexity_unimodal(:,h,handles.curlabel);
            %SEM = grd.best_Complexity_unimodal(:,h,handles.curlabel);
            ylb = 'Model complexity';

    end
end

% Ps              = handles.MLparams.ParamCombs;
% Ps_desc         = handles.MLparams.ParamDesc;
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
    cpt = {'g'};
end

x = repmat((1:lPi)',nPi,1); y = P(:); hi = zeros(nPi,1);
axes(handles.axes17); cla; hold on

mlx = Inf; mux = 0;
% Plot performance graphs
xi      = (1:lPi)'; 
hi      = bar(xi,P,'LineWidth', 1.0);
mlxi    = nm_nanmin(P,[],2); muxi = nm_nanmax(P,[],2);
for i=1:nPi, hi(i).FaceColor = cpt(i,:); end
for i=1:lPi, axl{i} = sprintf('M%g',i); end
% Plot markers
%plot(x, y, 'k.', 'MarkerSize', 1)
xlim(handles.axes17,'auto')
ylim(handles.axes17,'auto')
xticks(handles.axes17,xi);
xticklabels(handles.axes17,axl);

switch meastype
    case 'Cross-validation performance summary'
        hl = plot([xi-0.5 xi+0.5],repmat(grd.best_TSperf_multimodal(h),lPi,1),'k--','LineWidth',2);
    case 'Multi-class cross-validation performance summary'
        if strcmpi(pardesc.GridParam,'BAC')
            hl = plot([xi-0.5 xi+0.5],repmat(grd.MultiClass.BAC,lPi,1),'k--','LineWidth',2);
        else
            hl = plot([xi-0.5 xi+0.5],repmat(grd.MultiClass.accuracy,lPi,1),'k--','LineWidth',2);
        end
end 
if exist('lgstr','var') && ~isempty(lgstr) 
    handles.legend_modelperf = legend(lgstr,'Location','Best','FontSize',handles.LegendFontSize,'LineWidth',1); 
else
    legend('hide')
end
% Set GUI data
% pss = cell(1,numel(y));
% ll=1;
% for k=1:size(P,2)
%     for i=1:size(P,1)
%         psj=sprintf(' %s: %1.3f', ylb, P(i,k));
%         for j = 1:size(Ps_desc{h},2)
%             if iscell(Ps{h}(i,j)),
%                 ijPs = Ps{h}{i,j};
%             else
%                 ijPs = Ps{h}(i,j);
%             end
%             if ~ischar(ijPs), ijPs = num2str(ijPs); end
%             psj = sprintf('%s\n%s: %s',psj, Ps_desc{h}{j}, ijPs);
%         end
%         psj = psj(2:end);
%         pss{ll} = sprintf('%s\n',psj);
%         ll=ll+1;
%     end
% end
% figdata.x = x;
% figdata.y = y;
% figdata.patterntext = repmat(pss', nPi,1);
% figdata.parentui = handles.pnModelPerf;
% figdata.textHdl = annotation(figdata.parentui, 'textbox', 'String','', ...
%                             'Interpreter','none', ... %'VerticalAlign', 'Top', ...
%                             'Color', 'black', ...
%                             'BackgroundColor',[.7 .9 .7], ...
%                             'EdgeColor', 'black', ...
%                             'LineWidth', 0.5, ...
%                             'Margin', 5, ...
%                             'FitBoxToText','on', ...
%                             'Visible','off');
% set(handles.axes17,'UserData',figdata);

% Define x and y axis label
% strp = [];
% for i = 1 : numel(Ps_desc{h})
%     strp = sprintf('%s, %s', strp,Ps_desc{h}{i} ); 
% end
%strp = strp(3:end);
xlabel('Unimodal classifiers');
ylabel(ylb);
%xlim([0 size(P,1)]);
cla(handles.axes35); cla(handles.axes37);
handles.axes17.Position(3) = handles.axes37.Position(1)+ handles.axes37.Position(3) - handles.axes17.Position(1) ;
handles.axes35.Visible = 'off'; handles.axes35.YLabel.Visible = 'off'; handles.axes35.XLabel.Visible = 'off'; handles.axes35.Title.Visible = 'off'; 
handles.axes37.Visible = 'off'; handles.axes37.YLabel.Visible = 'off'; handles.axes37.XLabel.Visible = 'off'; handles.axes37.Title.Visible = 'off'; 
handles.cmdExportModelPerfC1.Visible='off'; handles.cmdExportModelPerfC2.Visible='off';
   
handles.cmdExportModelPerfVec.Position(1) = handles.axes17.Position(1)+handles.axes17.Position(3)+0.005;