function load_selModality(handles)

if strcmp(handles.popupmenu1.String{handles.popupmenu1.Value},'Multi-group classifier')
    multiflag = true;
else
    multiflag = false;
end

popuplist=[];
if ~multiflag

    for i=1:numel(handles.visdata{handles.curlabel})
        popuplist{i} = sprintf('Modality %g: %s', i, handles.NM.datadescriptor{handles.visdata{handles.curlabel}{handles.curmodal}.params.varind}.desc);
    end
    handles.selModality.String = popuplist;
    v = handles.visdata{handles.curlabel}{handles.curmodal};
    popuplist=[];
    if isfield(v,'MEAN'),                      popuplist{1} = 'Feature weights [Overall Mean (StErr)]';                        end
    if isfield(v,'MEAN_CV2'),                  popuplist{end+1} = 'Feature weights [Grand Mean (StErr)]';                      end
    if isfield(v,'CVRatio'),                   popuplist{end+1} = 'CV-ratio of feature weights [Overall Mean]';                end
    if isfield(v,'CVRatio_CV2'),               popuplist{end+1} = 'CV-ratio of feature weights [Grand Mean]';                  end
    if isfield(v,'FeatProb'),                  popuplist{end+1} = 'Feature selection probability [Overall Mean]';              end
    if isfield(v,'Prob_CV2'),                  popuplist{end+1} = 'Probability of feature reliability (95%-CI) [Grand Mean]';  end
	if isfield(v,'SignBased_CV2'),             popuplist{end+1} = 'Sign-based consistency';  									end
	if isfield(v,'SignBased_CV2_z'),           popuplist{end+1} = 'Sign-based consistency (Z score)';  						end
	if isfield(v,'SignBased_CV2_p_uncorr'),    popuplist{end+1} = 'Sign-based consistency -log10(P value)';  					end
	if isfield(v,'SignBased_CV2_p_fdr'),       popuplist{end+1} = 'Sign-based consistency -log10(P value, FDR)';  				end
    if isfield(v,'Spearman_CV2'),              popuplist{end+1} = 'Spearman correlation [Grand Mean]';                         end
    if isfield(v,'Pearson_CV2'),               popuplist{end+1} = 'Pearson correlation [Grand Mean]';                          end
    if isfield(v,'Spearman_CV2_p_uncorr'),     popuplist{end+1} = 'Spearman correlation -log10(P value) [Grand Mean]';         end
    if isfield(v,'Pearson_CV2_p_uncorr'),      popuplist{end+1} = 'Pearson correlation -log10(P value) [Grand Mean]';          end
    if isfield(v,'Spearman_CV2_p_fdr'),        popuplist{end+1} = 'Spearman correlation -log10(P value, FDR) [Grand Mean]';    end
    if isfield(v,'Pearson_CV2_p_fdr'),         popuplist{end+1} = 'Pearson correlation -log10(P value, FDR) [Grand Mean]';     end
    if isfield(v,'PermProb_CV2'),              popuplist{end+1} = 'Permutation-based -log10(P value) [Grand Mean]';            end
    if isfield(v,'PermProb_CV2_FDR_PVAL'),     popuplist{end+1} = 'Permutation-based -log10(P value, FDR) [Grand Mean]';       end
    if isfield(v,'PermZ_CV2'),                 popuplist{end+1} = 'Permutation-based Z Score [Grand Mean]';                    end
	if isfield(v,'Analytical_p'),     			popuplist{end+1} = 'Analytical -log10(P Value) for Linear SVM [Grand Mean]';    end
	if isfield(v,'Analyitcal_p_fdr'),     		popuplist{end+1} = 'Analytical -log10(P Value, FDR) for Linear SVM [Grand Mean]';end
    if isfield(v,'PermModel_Eval_Global'),     popuplist{end+1} = 'Model P value histogram';                                   end
    handles.selVisMeas.String = popuplist; 
    VisOnFl = 'on';
    VisElFl = 'on';
else
    if isfield(v,'ObsModel_Eval_Global_Multi')
        handles.selVisMeas.Enable = 'on'; 
        popuplist{1} = 'Model P value histogram [Multi-group]';
        for i=1:numel(v.ObsModel_Eval_Global_Multi_Bin)
            popuplist{end+1} = sprintf('Model P value histogram [Class %s vs. Rest]',handles.NM.groupnames{i});    
        end
        handles.selVisMeas.String = popuplist;
        VisOnFl = 'on';
        VisElFl = 'off';
    else
        VisOnFl = 'on';
        VisElFl = 'on';
    end
end
% Toggle enabled status of visualization controls
handles.selVisMeas.Enable = VisOnFl; 
handles.selModality.Enable = VisElFl;
handles.selPager.Enable = VisElFl;
handles.tglSortFeat.Enable = VisElFl;
handles.cmdExportFeats.Enable = VisElFl;