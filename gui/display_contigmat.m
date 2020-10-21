% =========================================================================
% =                   CONTINGENCY MATRIX INFO                             = 
% =========================================================================
function handles = display_contigmat(handles, contigmat)

contig=[];
GraphType = get(handles.selYaxis,'Value');

if ~exist('contigmat','var') || isempty(contigmat)
    
    switch handles.modeflag

        case 'regression'

            contig{end+1} = ['\bf Coefficient of determination [%]: \rm ' num2str(handles.Regr.R2(handles.curlabel),'%1.1f')];
            contig{end+1} = ['\bf Pearson r (CI-95%): \rm' num2str(handles.Regr.r(handles.curlabel),'%1.2f') ...
                                ' (' num2str(handles.Regr.r_95CI_low(handles.curlabel),'%1.2f') '-' num2str(handles.Regr.r_95CI_up(handles.curlabel),'%1.2f')  ')'];
            contig{end+1} = ['\bf P(T) value: \rm' num2str(handles.Regr.p(handles.curlabel),'%g') ' (' num2str(handles.Regr.t(handles.curlabel),'%1.2f') ')' ];
            contig{end+1} = ['\bf Mean Absolute Error: \rm' num2str(handles.Regr.MAE(handles.curlabel),'%1.1f')];
            contig{end+1} = ['\bf Mean Squared Error: \rm' num2str(handles.Regr.MSE(handles.curlabel),'%1.1f')];
            contig{end+1} = ['\bf NRMSD: \rm' num2str(handles.Regr.NRSMD(handles.curlabel),'%1.1f')];
            contigmat = handles.Regr.contigmat;

        case 'classification'
            
            h_class  = get(handles.popupmenu1,'Value');
            h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
            h_classlist     = get(handles.popupmenu1,'String');
            if strcmpi(h_classlist{h_class},'Multi-group classifier') && h_onevsall_val > 1
                contigmat = handles.MultiClass.class{h_onevsall_val-1};
            else
                switch GraphType
                    case {4,5,6}
                         contigmat = handles.BinClass{h_class}.prob_contingency;
                    otherwise
                        contigmat = handles.BinClass{h_class}.contingency;
                end
            end
    end
end
contig{end+1} = ['\bf TP / TN:\rm ' num2str(contigmat.TP,'%1.0f') ' / ' num2str(contigmat.TN,'%1.0f')];
contig{end+1} = ['\bf FP / FN:\rm ' num2str(contigmat.FP,'%1.0f') ' / ' num2str(contigmat.FN,'%1.0f')];
contig{end+1} = ['\bf Accuracy [%]:\rm ' num2str(contigmat.acc,'%1.1f')];
contig{end+1} = ['\bf Sensitivity [%]:\rm ' num2str(contigmat.sens,'%1.1f')];
contig{end+1} = ['\bf Specificity [%]:\rm ' num2str(contigmat.spec,'%1.1f')];
contig{end+1} = ['\bf Balanced Accuracy [%]:\rm ' num2str(contigmat.BAC,'%1.1f')];
if isfield(contigmat,'AUC')
    contig{end+1} = ['\bf Area Under the Curve:\rm ' num2str(contigmat.AUC,'%1.2f')];
end
contig{end+1} = ['\bf Matthews Coorelation Coefficient:\rm ' num2str(contigmat.MCC,'%1.1f')];
contig{end+1} = ['\bf Positive Predictive Value [%]:\rm ' num2str(contigmat.PPV,'%1.1f')];
contig{end+1} = ['\bf Negative Predictive Value [%]:\rm ' num2str(contigmat.NPV,'%1.1f')];
contig{end+1} = ['\bf False Positive Rate:\rm ' num2str(contigmat.FPR,'%1.1f')];
contig{end+1} = ['\bf Positive Likelihood Ratio:\rm ' num2str(contigmat.pLR,'%1.1f')];
contig{end+1} = ['\bf Negative Likelihood Ratio:\rm ' num2str(contigmat.nLR,'%1.1f')];
contig{end+1} = ['\bf Prognostic Summary Index:\rm ' num2str(contigmat.PSI,'%1.1f')];
contig{end+1} = ['\bf Youden''s J statistic:\rm ' num2str(contigmat.Youden,'%1.1f')];
contig{end+1} = ['\bf # Needed to Predict/Diagnose:\rm ' num2str(contigmat.NNP,'%1.1f') '/' num2str(contigmat.NND,'%1.1f')];
contig{end+1} = ['\bf Diagnostic Odds Ratio:\rm ' num2str(contigmat.DOR,'%1.1f')];

cla(handles.axes5);
delete(findall(handles.figure1,'Tag','AnnotPerfMeas'))

switch handles.NM.modeflag
    case 'classification'
        y_start = 0.5;
        y_height = 0.45;
    case 'regression'
        y_start = 0.39;
        y_height = 0.58;
end
handles.txtPerf = annotation(handles.pnContigCmds, ...
                        'textbox',[handles.axes20.Position(1) y_start handles.axes20.Position(3) y_height] , ...
                        'String',contig, ...
                        'FitBoxToText','off', ...
                        'FontUnits', 'normalized', ...
                        'FontSize', 0.0185, ...
                        'Margin', 5, ...
                        'Units', 'normalized', ...
                        'LineWidth', 1.5, ...
                        'Tag','AnnotPerfMeas', ...
                        'Interpreter','tex');

end