function OOCVinfo = nk_GetOOCVInfo(NM, searchmode)

if ~exist('searchmode','var') || isempty(searchmode), searchmode = 'oocv'; end

OOCVinfo.searchmode = searchmode;
AS = nk_GetAnalysisStatus(NM);

if AS.analcompleteflag >0 && AS.oocvappflag
    
    OOCVinfo.AnalVec = find(AS.completed_analyses);
    OOCVinfo.n_AnalCompl = numel(OOCVinfo.AnalVec);
    OOCVinfo.n_OOCV = numel(NM.OOCV);
    
    switch searchmode

        case 'oocv'

            for i = 1 : OOCVinfo.n_OOCV
                OOCVinfo.desc{i} = NM.OOCV{i}.desc;
                OOCVinfo.date{i} = NM.OOCV{i}.date;
                OOCVinfo.n_subjects_all(i) = NM.OOCV{i}.n_subjects_all;
                OOCVinfo.labels_known(i) = NM.OOCV{i}.labels_known;
                if OOCVinfo.labels_known(i)
                    OOCVinfo.n_subjects{i} = NM.OOCV{i}.n_subjects;
                    OOCVinfo.label{i} = NM.OOCV{i}.label;
                else
                    OOCVinfo.n_subjects{i} = [];
                    OOCVinfo.label{i} = [];
                end
                OOCVinfo.cases{i} = NM.OOCV{i}.cases;
                OOCVinfo.files{i} = NM.OOCV{i}.files;
                OOCVinfo.defs{i} = NM.OOCV{i}.defs;
                
                for j = 1: OOCVinfo.n_AnalCompl 
                    OOCVinfo.Analyses{j}.Num = OOCVinfo.AnalVec(j);
                    if isfield(NM.analysis{OOCVinfo.AnalVec(j)},'OOCV') && ~isempty(NM.analysis{OOCVinfo.AnalVec(j)}.OOCV{i})
                        OOCVinfo.Analyses{j}.OOCVavail = true;
                        OOCVinfo.Analyses{j}.predictions = NM.analysis{OOCVinfo.AnalVec(j)}.OOCV{i}.predictions;
                        if isfield(NM.analysis{OOCVinfo.AnalVec(j)}.OOCV{i},'multi_predictions')
                            OOCVinfo.Analyses{j}.multi_predictions = NM.analysis{OOCVinfo.AnalVec(j)}.OOCV{i}.multi_predictions;
                        end
                    else
                        OOCVinfo.Analyses{j}.OOCVavail = false;
                    end
                end
            end

        case 'analysis'

             for i = 1 : OOCVinfo.n_AnalCompl
                 OOCVinfo.Analyses{i}.num = OOCVinfo.AnalVec(i);
                 if isfield(NM.analysis{OOCVinfo.AnalVec(i)},'OOCV')
                    OOCVinfo.Analyses{i}.OOCVdone = true;
                    OOCVinfo.Analyses{i}.OOCVvec = find(~cellfun(@isempty,NM.analysis{OOCVinfo.AnalVec(i)}.OOCV));
                    for j = 1:numel(OOCVinfo.Analyses{i}.OOCVvec)
                        try
                            OOCVinfo.Analyses{i}.desc{j} = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.desc;
                            OOCVinfo.Analyses{i}.date{j} = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.date;
                            OOCVinfo.Analyses{i}.n_subjects_all(j) = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.n_subjects_all;
                            OOCVinfo.Analyses{i}.labels_known(j) = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.labels_known;

                            if OOCVinfo.Analyses{i}.labels_known(j)
                                lblstr = sprintf(', labels known');
                                OOCVinfo.Analyses{i}.n_subjects{j} = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.n_subjects;
                                OOCVinfo.Analyses{i}.label{j} = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.label;
                            else
                                OOCVinfo.Analyses{i}.n_subjects{j} = [];
                                lblstr = '';
                                 OOCVinfo.Analyses{i}.label = [];
                            end
                            OOCVinfo.Analyses{i}.cases{j}   = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.cases;
                            OOCVinfo.Analyses{i}.files{j}   = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.files;
                            OOCVinfo.Analyses{i}.defs{j}    = NM.OOCV{OOCVinfo.Analyses{i}.OOCVvec(j)}.defs;
                            OOCVinfo.Analyses{i}.descriptor{j} = sprintf('%s (<-%s): %g cases%s',  OOCVinfo.Analyses{i}.desc{j}, OOCVinfo.Analyses{i}.date{j}, OOCVinfo.Analyses{i}.n_subjects_all(j), lblstr); 
                        catch
                            cprintf('red','\nOOCV data container %g does not exist',j);
                        end
                    end
                 else
                     OOCVinfo.Analyses{i}.OOCVdone = false;
                 end
             end

    end
    
end