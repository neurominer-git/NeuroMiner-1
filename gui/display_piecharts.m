% =========================================================================
% =                          PIE CHARTS                                   =
% =========================================================================
function [h1pie, h2pie] = display_piecharts(handles, Pre1, contigmat, labelh)
GraphType = get(handles.selYaxis,'Value');
% Prepare the pie chart descriptors

switch handles.modeflag
    case 'classification'
        h_class = get(handles.popupmenu1,'Value');
        h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
        h_classlist     = get(handles.popupmenu1,'String');
        if strcmpi(h_classlist{h_class},'Multi-group classifier') 
            if h_onevsall_val > 1
                colind = h_onevsall_val-1;
                groupnames = {handles.NM.groupnames{colind}, 'REST'};
                if ~exist('contigmat','var') || isempty(contigmat)
                    contigmat = handles.MultiClass.class{colind};
                    labelh = handles.MultiClass.onevsall_labels(:,colind);
                end
                cl1 = handles.colptin(colind,:);
                cl2 = [0.6 0.6 0.6];
            end
        else
            groupnames = handles.BinClass{h_class}.groupnames;
            if ~exist('contigmat','var') || isempty(contigmat)
                switch GraphType
                    case {4,5,6}
                         contigmat = handles.BinClass{h_class}.prob_contingency;
                    otherwise
                        contigmat = handles.BinClass{h_class}.contingency;
                end
                labelh = handles.BinClass{h_class}.labelh;
            end
            cl1 = handles.colptin(handles.BinClass{h_class}.groupind(1),:);
           
            if ~handles.BinClass{h_class}.one_vs_all
                cl2 = handles.colptin(handles.BinClass{h_class}.groupind(2),:);
            else
                cl2= rgb('DarkGrey');
            end
        end
    case 'regression'
        set(handles.txtPretestProb,'visible','on');
        thresh = get(handles.txtBinarize,'String');
        groupnames{1} = sprintf('Group 1 (>=%s)',thresh);
        groupnames{2} = sprintf('Group 2 (< %s)',thresh);
        if ~exist('contigmat','var') || isempty(contigmat)
            contigmat = handles.Regr.contigmat;
            labelh = handles.Regr.b_label;
        end
        cl1 = [0.5 0.5 0.5]; cl2 = [0.5 0.5 0.5];
end
if ~exist('Pre1','var') || isempty(Pre1),
    Pre1    = sum(labelh==1)/numel(labelh);
end

% Display pie charts
LR1     = contigmat.pLR;
h1pie = display_pie(handles.axes3, handles.cmdExportPies, Pre1, LR1 , groupnames{1}, cl1);

Pre2    = 1-Pre1;
LR2     = (contigmat.spec/100) / (1-(contigmat.sens/100));
h2pie = display_pie(handles.axes4, handles.cmdExportPies, Pre2, LR2 , groupnames{2}, cl2);

if isempty(h1pie) || isempty(h2pie)
    set(handles.txtPretestProb,'Visible','off')
else
    set(handles.txtPretestProb,'Visible','on')
end

function hpie = display_pie(axesh, cmdh, Pre, LR , groupnames, cl)

Post = Pre*LR/(1-Pre+Pre*LR); 
if ~isfinite(Post),Post=1; end
Gain = Post-Pre;

axes(axesh); cla; 
if ~isnan(Post)
    axesh.Visible = 'on';cmdh.Visible = 'on';
    if Gain>0
        hpie   = pie(axesh,[Pre Gain],{sprintf('Pre\n%1.1f%%',Pre*100),sprintf('Gain\n%1.1f%%',(Post-Pre)*100)});   
    else
        hpie   = pie(axesh,Pre,{sprintf('Pre\n%1.1f%%',Pre*100)});
    end
    hp      = findobj(hpie, 'Type', 'patch');
    set(hp(1), 'FaceColor', cl, 'FaceAlpha', 0.3);
    if numel(hp)>1; set(hp(2), 'FaceColor', cl, 'FaceAlpha', 1); end
    title(groupnames); axesh.Title.Position(2) = 1.1; axesh.FontUnits = 'normalized'; axesh.FontSize=0.1;
    display_patch_labels(hpie, Pre, Gain)
else
    axesh.Visible='off'; cmdh.Visible = 'off';
    hpie=[];
end

function display_patch_labels(pie, pre, gain)

textPos1 = pie(2).Position; 
pie(2).FontWeight = 'bold';
if numel(pie)>2; pie(4).FontWeight = 'bold'; end
if pre>0.15
    textPos1 = textPos1*0.1; % add offset
    pie(2).Color = 'k'; 
else
    %textPos1 = textPos1*1.05; % add offset
    pie(2).Color = 'k';
end
pie(2).Position = textPos1;
if gain>0 && numel(pie)>2
    textPos2 = pie(4).Position;
    if gain>0.15
        textPos2 = textPos2*0.1; % add offset
        pie(4).Color = 'w';
    elseif gain>0        
        %textPos2 = textPos2*1.05; % add offset
        pie(4).Color = 'k';
    end
    pie(4).Position = textPos2; 
end

    