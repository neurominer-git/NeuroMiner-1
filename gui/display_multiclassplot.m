% =========================================================================
% =                  MULTI-GROUP CLASSIFICATION PLOTS                     =
% =========================================================================
function handles = display_multiclassplot(handles)

axes(handles.axes1);
%uistack(handles.axes1,'top')
cla;
hold on

h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
probabilities   = handles.MultiClass.probabilities;
errors          = logical(handles.MultiClass.errors');
[m,n]           = size(probabilities);
subjvec         = (1:m)';
errx            = find(errors == 1);
erry            = zeros(size(errx));
[~,maxI]        = max(probabilities,[],2);
for i=1:length(errx)
    erry(i) = probabilities(errx(i),maxI(errx(i)));
end
ind0 = false(m,n); ind0(sub2ind([m n], 1:m, maxI')) = true;

for h=1:n
    indh = maxI == h; %subjvech = subjvec(indh); %ind0h = ind0(indh,:);
    if ~isempty(indh)
        
        if h_onevsall_val >1 && h ~= h_onevsall_val-1
            hcolptin = [0.7 0.7 0.7];
            %hcolpt = hcolptin;
            %hcolpx = hcolpt;
        else
            hcolptin = handles.colptin(h,:);
            %hcolpt = handles.colpt{h};
            %hcolpx = handles.colpx{h};
        end
        
        handles.MultiClass.bx(h)   = plot(subjvec,probabilities(:,h), ...
                                        handles.colpt, ...
                                        'MarkerSize',4, ...
                                        'MarkerFaceColor',hcolptin,...
                                        'Color',hcolptin, ...
                                        'LineWidth',1, ...
                                        'LineStyle','none');
        indhx =  indh & ~errors; indh0 = indh & errors;
        if sum(indhx)>0
            handles.MultiClass.b(h)    = plot(subjvec(indhx),probabilities(indhx,h), ...
                                        handles.colpt, ...
                                        'MarkerSize',handles.DataMarkerSize,...
                                        'MarkerFaceColor',hcolptin,...
                                        'LineWidth',handles.DataMarkerWidth,...
                                        'Color',hcolptin, ...
                                        'LineStyle','none');
        else
            handles.MultiClass.b(h) =plot(0,0);
        end
        if sum(indh0)>0
            handles.MultiClass.be(h)    = plot(subjvec( indh0 ),probabilities(indh0,h), ...
                                        handles.colpx, ...
                                        'MarkerSize',14,...
                                        'MarkerFaceColor',hcolptin,...
                                        'Color',hcolptin, ...
                                        'LineWidth',handles.DataMissMarkerWidth,...
                                        'LineStyle','none');
        else
            handles.MultiClass.be(h) = plot(0,0);
        end
        groupnames_err{h}           = 'misclassified';
    else
        handles.MultiClass.b(h)    = plot (0,0);
        handles.MultiClass.bx(h)    = plot (0,0);
    end
end

xlim([min(subjvec) max(subjvec)]);
%handles.MultiClass.be = plot(errx,erry,'kx','MarkerSize',handles.DataMissMarkerSize,'LineWidth',handles.DataMissMarkerWidth);
legend([handles.MultiClass.b handles.MultiClass.be], [handles.NM.groupnames(:)', groupnames_err(:)'], 'Location','Best', 'FontSize',handles.LegendFontSize)

% Define textbox info data 
pss = cell(1,m); subjects = handles.MultiClass.cases;
figdata.y = zeros(m,1);
for i=1:numel(pss)
    expgroupi   = handles.NM.groupnames{handles.MultiClass.labels(i)}; 
    predgroupi  = handles.NM.groupnames{maxI(i)}; 
    pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s'], subjects{i}, expgroupi, predgroupi);
    for j=1:handles.ngroups
        pss{i} = sprintf('%s\nProbability %s: %1.2f', pss{i}, handles.NM.groupnames{j}, probabilities(i,j));
    end
    figdata.y(i) = probabilities(i,maxI(i));
end
hText = uicontrol('Style','text','String', pss{numel(pss)},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.x           = subjvec;
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
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);
ylabel('Multi-class probabilities'); 
