% =========================================================================
% =                        VISUALIZATION PLOT                             =
% =========================================================================
function handles = display_visual(handles)
global st
st.ParentFig = handles.figure1;
varind = get(handles.selModality,'Value');
measind = get(handles.selVisMeas,'Value');
meas = get(handles.selVisMeas,'String');
load_selPager(handles)
pageind = get(handles.selPager,'Value');
page = get(handles.selPager,'String');
sortfl = get(handles.tglSortFeat,'Value');
axes(handles.axes33); cla; hold on
set(handles.axes33,'TickLabelInterpreter','none')

v = handles.visdata{handles.curlabel}{varind};
if v.params.visflag == 1
    featind = 1:v.params.nfeats;
else
    try
        featind = eval(page{pageind});
    catch
        featind = 1:v.params.nfeats;
    end
end

x = featind(1) - 0.5: featind(end);
curclass = get(handles.popupmenu1,'Value');

if strcmp(handles.popupmenu1.String{curclass},'Multi-group classifier')
    multiflag = true;
else
    multiflag = false;
end

if measind>numel(meas)
    measind=numel(meas);
    handles.selVisMeas.Value=measind;
end

switch meas{measind} 
    case 'Model P value histogram'
        fl = 'off';
    otherwise
        if isfield(v,'PermModel_Crit_Global_Multi') && multiflag
            fl = 'off';
        else
            fl = 'on';
        end
end

handles.selPager.Enable = fl;
handles.tglSortFeat.Enable = fl;
handles.cmdExportFeats.Enable = fl;
vlineval = [];
switch meas{measind}
    
    case {'Feature weights [Overall Mean (StErr)]','Feature weights [Grand Mean (StErr)]'}
        switch meas{measind}
            case 'Feature weights [Overall Mean (StErr)]'
                if iscell(v.MEAN)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.MEAN,[],2),2);
                        e = nm_nanmean(nk_cellcat(v.SE,[],2),2);
                    else
                        y = v.MEAN{curclass}; e = v.SE{curclass}; 
                    end
                else
                    y = v.MEAN; e = v.SE; 
                end
            case 'Feature weights [Grand Mean (StErr)]'
                if iscell(v.MEAN_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.MEAN_CV2,[],2),2);
                        e = nm_nanmean(nk_cellcat(v.SE_CV2,[],2),2);
                    else
                        y = v.MEAN_CV2{curclass}; e = v.SE_CV2{curclass}; 
                    end
                else
                    y = v.MEAN_CV2; e = v.SE_CV2; 
                end
        end
        y(~isfinite(y))=0;
        vals = y + (sign(y) .* e);
        miny = nanmin(vals);  maxy = nanmax(vals);
        if sortfl,
            [~, ind] = sort(abs(y),'descend');
            y = y(ind); y = y(featind);
            se = e(ind); se = se(featind);
        else
            ind = (1:v.params.nfeats)';
            y = y(featind);
            se = e(featind);
        end 
        
        switch v.params.visflag
            case {0, 3, 4, 5, 'matrix'}
                set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on');
                barh(x, y,'FaceColor','r');
                h = herrorbar(y, x , se, se,'ko');
                set(h,'MarkerSize',0.001);
            case 1
                st.fig = handles.pn3DView;
                st.NMaxes = [ handles.axes26 handles.axes27 handles.axes28]; 
                set(handles.pn3DView,'Visible','on'); set(handles.axes33,'Visible','off');
                nk_WriteVol(y,'temp',2,v.params.brainmask,v.params.badcoords,0,'gt');
                if ~isfield(handles,'orthviews'), 
                    handles.orthviews = nk_orthviews('Image','temp.nii'); 
                    colormap(jet);
                else
                    nk_orthviews('Redraw')
                end
        end
        
    otherwise
        switch meas{measind}
            
            case 'CV-ratio of feature weights [Overall Mean]'
                if iscell(v.CVRatio)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.CVRatio,[],2),2);
                    else
                        y = v.CVRatio{curclass}; 
                    end
                else
                    y = v.CVRatio; 
                end
                miny = nanmin(y); maxy = nanmax(y);
                
            case 'CV-ratio of feature weights [Grand Mean]'
                if iscell(v.CVRatio_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.CVRatio_CV2,[],2),2);
                    else
                        y = v.CVRatio_CV2{curclass};  
                    end
                else
                    y = v.CVRatio_CV2;
                end
                miny = nanmin(y); maxy = nanmax(y);
                
            case 'Feature selection probability [Overall Mean]'
                if iscell(v.FeatProb)
                    if multiflag
                        y = nm_nanmean(v.FeatProb{1},2);
                    else
                        y = v.FeatProb{1}(:,curclass);  
                    end   
                else
                    y = v.FeatProb;
                end
                miny = 0; maxy = nanmax(y);
                
            case 'Probability of feature reliability (95%-CI) [Grand Mean]'
                if iscell(v.Prob_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Prob_CV2,[],2),2);
                    else
                        y = v.Prob_CV2{curclass};  
                    end
                else
                    y = v.Prob_CV2;
                end
                miny = -1; maxy = 1;
                
            case 'Sign-based consistency'
                if iscell(v.SignBased_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.SignedBased_CV2,[],2),2);
                    else
                        y = v.SignBased_CV2{curclass};
                    end
                else
                    y = v.SignBased_CV2;
                end
                miny = 0; maxy = 1;
                
            case 'Sign-based consistency (Z score)'
                if iscell(v.SignBased_CV2_z)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.SignBased_CV2_z,[],2),2);
                    else
                        y = v.SignBased_CV2_z{curclass};
                    end
                else
                    y = v.SignBased_CV2_z;
                end
                miny = nanmin(y(:)); maxy = nanmax(y(:));
                
            case 'Sign-based consistency -log10(P value)'
                if iscell(v.SignBased_CV2_p_uncorr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.SignBased_CV2_p_uncorr,[],2),2);
                    else
                        y = v.SignBased_CV2_p_uncorr{curclass};
                    end
                else
                    y = v.SignBased_CV2_p_uncorr;
                end 
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
            case 'Sign-based consistency -log10(P value, FDR)'
                if iscell(v.SignBased_CV2_p_fdr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.SignBased_CV2_p_fdr,[],2),2);
                    else
                        y = v.SignBased_CV2_p_fdr{curclass};
                    end
                else
                    y = v.SignBased_CV2_p_fdr;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
            case 'Spearman correlation [Grand Mean]'
                if iscell(v.Spearman_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Spearman_CV2,[],2),2);
                    else
                        y = v.Spearman_CV2{curclass};
                    end
                else
                    y = v.Spearman_CV2;
                end
                miny = nanmin(y(:)); maxy = nanmax(y(:));
            case 'Pearson correlation [Grand Mean]'
                if iscell(v.Pearson_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Pearson_CV2,[],2),2);
                    else
                        y = v.Pearson_CV2{curclass};  
                    end
                else
                    y = v.Pearson_CV2;
                end
                miny = nanmin(y(:)); maxy = nanmax(y(:));
                
            case 'Spearman correlation -log10(P value) [Grand Mean]'
                if iscell(v.Spearman_CV2_p_uncorr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Spearman_CV2_p_uncorr,[],2),2);
                        if isfield (v,'Spearman_CV2_p_uncorr_STD')
                            se = nm_nanmean(sqrt(nk_cellcat(v.Spearman_CV2_p_uncorr_STD,[],2)),2);
                        end
                    else
                        y = v.Spearman_CV2_p_uncorr{curclass}; 
                        if isfield (v,'Spearman_CV2_p_uncorr_STD')
                            se = sqrt(v.Spearman_CV2_p_uncorr_STD{curclass});
                        end
                    end
                else
                    y = v.Spearman_CV2_p_uncorr;
                    if isfield (v,'Spearman_CV2_p_uncorr_STD')
                        se = sqrt(v.Spearman_CV2_p_uncorr_STD);
                    end
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Pearson correlation -log10(P value) [Grand Mean]'
                if iscell(v.Pearson_CV2_p_uncorr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Pearson_CV2_p_uncorr,[],2),2);
                        if isfield (v,'Pearson_CV2_p_uncorr_STD')
                            se = nm_nanmean(sqrt(nk_cellcat(v.Pearson_CV2_p_uncorr_STD,[],2)),2);
                        end
                    else
                        y = v.Pearson_CV2_p_uncorr{curclass};  
                        if isfield (v,'Pearson_CV2_p_uncorr_STD')
                            se = sqrt(v.Pearson_CV2_p_uncorr_STD{curclass});
                        end
                    end
                else
                    y = v.Pearson_CV2_p_uncorr;
                    if isfield (v,'Pearson_CV2_p_uncorr_STD')
                        se = sqrt(v.Pearson_CV2_p_uncorr_STD);
                    end
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Spearman correlation -log10(P value, FDR) [Grand Mean]'
                if iscell(v.Spearman_CV2_p_fdr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Spearman_CV2_p_fdr,[],2),2);
                    else
                        y = v.Spearman_CV2_p_fdr{curclass}; 
                    end
                else
                    y = v.Spearman_CV2_p_fdr;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Pearson correlation -log10(P value, FDR) [Grand Mean]'
                if iscell(v.Pearson_CV2_p_fdr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Pearson_CV2_p_fdr,[],2),2);
                    else
                        y = v.Pearson_CV2_p_fdr{curclass};  
                    end
                else
                    y = v.Pearson_CV2_p_fdr;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Permutation-based Z Score [Grand Mean]'
                if iscell(v.PermZ_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.PermZ_CV2,[],2),2);
                    else
                        y = v.PermZ_CV2{curclass};  
                    end
                else
                    y = v.PermZ_CV2;
                end
                miny = nanmin(y(:)); maxy = nanmax(y(:));
                
            case 'Permutation-based -log10(P value) [Grand Mean]'
                if iscell(v.PermProb_CV2)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.PermProb_CV2,[],2),2);
                    else
                        y = v.PermProb_CV2{curclass};  
                    end
                else
                    y = v.PermProb_CV2;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Permutation-based -log10(P value, FDR) [Grand Mean]'
                if iscell(v.PermProb_CV2_FDR)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.PermProb_CV2_FDR,[],2),2);
                    else
                        y = v.PermProb_CV2_FDR{curclass};  
                    end
                else
                    y = v.PermProb_CV2_FDR;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Analytical -log10(P Value) for Linear SVM [Grand Mean]'
                if iscell(v.Analytical_p)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Analytical_p,[],2),2);
                    else
                        y = v.Analytical_p{curclass};  
                    end
                else
                    y = v.Analytical_p;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Analytical -log10(P Value, FDR) for Linear SVM [Grand Mean]'
                if iscell(v.Analyitcal_p_fdr)
                    if multiflag
                        y = nm_nanmean(nk_cellcat(v.Analyitcal_p_fdr,[],2),2);
                    else
                        y = v.Analyitcal_p_fdr{curclass};  
                    end
                else
                    y = v.Analyitcal_p_fdr;
                end
                miny = 0; maxy = nanmax(y(:));
                vlineval = -log10(0.05);
                
            case 'Model P value histogram'
                if multiflag
                    y = nm_nanmean(v.PermModel_Crit_Global); 
                    vp = nm_nanmean(v.ObsModel_Eval_Global);
                    ve = nm_nanmean(v.PermModel_Eval_Global);
                else
                    y = v.PermModel_Crit_Global(curclass,:); 
                    vp = v.ObsModel_Eval_Global(curclass);
                    ve = v.PermModel_Eval_Global(curclass,:);
                end
                perms = length(v.PermModel_Eval_Global(1,:));
        end
        
        if multiflag && isfield(v,'PermModel_Crit_Global_Multi')
            if measind == 1
                 y = v.PermModel_Crit_Global_Multi; 
                 vp = v.ObsModel_Eval_Global_Multi;
                 ve = v.PermModel_Eval_Global_Multi;
            else
                 y = v.PermModel_Crit_Global_Multi_Bin(measind-1,:); 
                 vp = v.ObsModel_Eval_Global_Multi_Bin(measind-1);
                 ve = v.PermModel_Eval_Global_Multi_Bin(measind-1,:);
            end
            perms = length(v.PermModel_Eval_Global_Multi);
            meas{measind} = 'Model P value histogram';
        end
        
        y(~isfinite(y))=0;
        
        if ~strcmp(meas{measind},'Model P value histogram')
            if sortfl, 
                [~,ind] = sort(abs(y),'descend');
                y = y(ind); y = y(featind); 
                if exist('se','var'), se = se(ind); se=se(featind); end
            else
                ind = (1:v.params.nfeats)';
                y = y(featind);
                if exist('se','var'), se=se(featind); end
            end
            switch v.params.visflag
                case {0, 3, 4, 5}
                    set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on');
                    barh(x, y,'FaceColor','b','BarWidth',0.9,'FaceColor',rgb('SkyBlue'),'EdgeColor','w','LineWidth',1)
                    if exist('se','var')
                        h = herrorbar(y, x , se, se,'ko');
                        set(h,'MarkerSize',0.001);
                        ylimoff = nanmax(se);
                    else 
                        ylimoff = 0;
                    end
                    if y <= 0
                        ylim(handles.axes33,[nanmin(y)+ylimoff 0]);
                    else
                        ylim(handles.axes33,[0 nanmax(y)+ylimoff]);
                    end
                    if ~isempty(vlineval), xline(vlineval,'LineWidth',1.5, 'Color', 'red'); end
                case 1
                    st.fig = handles.pn3DView; 
                    st.NMaxes = [ handles.axes26 handles.axes27 handles.axes28];
                    set(handles.pn3DView,'Visible','on'); set(handles.axes33,'Visible','off');
                    nk_WriteVol(y,'temp',2,v.params.brainmask,[],0,'gt');
                    if ~isfield(handles,'orthviews'), 
                        handles.orthviews = nk_orthviews('Image','temp.nii'); 
                        %colormap(jet);
                    else
                        nk_orthviews('Redraw')
                    end

            end
            legend('off');
        else
             set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on');
             ah=histogram(handles.axes33,y,'Normalization','probability','BinWidth',2.5,'EdgeColor','none','FaceColor',rgb('CornflowerBlue')); 
             maxah= nanmax(ah.Values); ylim([0 maxah]); 
             handles.axes33.YTick = 0:maxah/10:maxah;
             yticklabels(handles.axes33,'auto')
             [xl,xlb]=nk_GetScaleYAxisLabel(handles.NM.analysis{handles.curranal}.params.TrainParam.SVM);
             xlim(xl); miny = xl(1); maxy=xl(2);
             xlabel(['Optimization criterion: ' xlb]);
             ylabel('Probability');
             hold on;
             xp = [ vp vp ]; yp = [ 0 maxah ];
             hl=line(xp, yp ,'LineWidth',2,'Color','r');
             Pval = sum(ve)/size(ve,2);
             if Pval ~= 0
                Pvalstr = sprintf('P=%g',Pval);
             else
                Pvalstr = sprintf('P<%g',1/perms);
             end
             legend(hl,Pvalstr);
        end
end

if  ~strcmp(meas{measind},'Model P value histogram')
    switch v.params.visflag
        case {0,3,4,5,'matrix'}
            xlabel(meas{measind},'FontWeight','bold');
            ylabel('Features','FontWeight','bold');
            set(gca,'YTick',x);
            if sortfl, 
                feats = v.params.features(ind); 
            else
                feats = v.params.features;
            end 
            feats = feats(featind);
            set(gca,'YTickLabel',feats);
            ylim([x(1)-0.5 x(end)+0.5]);
            if miny>=0,
                xlim([0 maxy]);
            else
                xlim([miny maxy]);
            end
           
            handles.visdata_table(handles.curlabel, handles.curmodal) = create_visdata_tables(v, handles.visdata_table(handles.curlabel, handles.curmodal), ind, 'reorder');
    end
end
guidata(handles.figure1, handles);