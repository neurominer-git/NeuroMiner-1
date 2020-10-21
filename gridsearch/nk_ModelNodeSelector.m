function [GD, MultiBinBind] = nk_ModelNodeSelector(GD, MD, label, f, d, nclass, Ps, Pdesc, combcell, act)

global CV MULTI GRD SAV MODEFL RFE RAND MULTILABEL SVM

MultiBinBind = [];

[op, ~, ~, PerfCritBin] = nk_ReturnEvalOperator(SVM.GridParam);
PerfCritMult=PerfCritBin;

for curlabel=1:MULTILABEL.dim

    %%%%%%%%%%%%%%%% SELECT PARAMS AT OPTIMUM %%%%%%%%%%%%%%%%%
    %% BINARY CLASSIFIERS / REGRESSORS

    % Optionally, apply penalization term on the grid search to avoid
    % overfitting and promote sparse and diverse models
    if isfield(GRD,'OptRegul') && GRD.OptRegul.flag && isfield(GRD.OptRegul,'type') 
        switch GRD.OptRegul.type
            case 1 % Model complexity is regularization constraint
                Regul = nk_ComputeRegulFunction(GD.C(:,:,curlabel), 0, GRD.OptRegul.RegulTypeComplexity);
            case 2
                % Model diversity
                if RFE.Filter.SubSpaceFlag && RFE.Filter.SubSpaceStrategy > 1
                    Regul = GD.M_DivV(:,:,curlabel);
                    for curclass = 1 : nclass
                        Regul(:,curclass) = nk_ComputeRegulFunction(GD.M_DivV(:,curclass,curlabel), 0, GRD.OptRegul.RegulTypeDiversity);
                    end
                else
                    Regul = nk_ComputeRegulFunction(GD.CV2Div(:,:,curlabel), 0, GRD.OptRegul.RegulTypeDiversity);
                end
            case 3
                % Combine Model diversity & complexity
                if RFE.Filter.SubSpaceFlag  && RFE.Filter.SubSpaceStrategy > 1
                    RegulDiv = GD.M_DivV(:,:,curlabel); RegulC = GD.C(:,:,curlabel);
                    for curclass = 1 : nclass
                        RegulDiv(:,curclass) = nk_ComputeRegulFunction(GD.M_DivV(:,curclass,curlabel), 0, GRD.OptRegul.RegulTypeDiversity);
                    end
                    for curclass = 1 : nclass
                        RegulC(:,curclass) = nk_ComputeRegulFunction(GD.C(:,curclass,curlabel), 0, GRD.OptRegul.RegulTypeComplexity);
                    end
                else
                    RegulDiv = nk_ComputeRegulFunction(GD.CV2Div(:,:,curlabel), 0, GRD.OptRegul.RegulTypeDiversity);
                    RegulC = nk_ComputeRegulFunction(GD.C(:,:,curlabel), 0, GRD.OptRegul.RegulTypeComplexity);
                end
                Regul = (RegulDiv + RegulC) / 2;
        end
    else
        Regul =  GD.C(:,:,curlabel);
    end

    % Loop through dichotomizers
    for curclass=1:nclass

        switch MODEFL

            case 'classification'
                if RAND.Decompose ~= 9
                    fprintf('\n\n');
                    cprintf('*black','**** Selected binary classifiers #%g (%s) ****', ...
                        curclass, CV.class{f,d}{curclass}.groupdesc)
                else
                    fprintf('\n\n');
                    cprintf('*black','**** Selected multi-group classifier ****')
                end
            case 'regression'
                fprintf('\n\n');
                cprintf('*black','**** Selected regression models ****')
        end

        if ~isfield(GD,'Weights'), GD.Weights = cell(size(GD.FEAT)); end
        % Now select optimal parameters for the learning machine
         nPerc = 1; PercVec = 0;
        if isfield(GRD,'NodeSelect') 
            switch GRD.NodeSelect.mode 
                case {2,3}
                    PercVec = GRD.NodeSelect.perc; nPerc = 1;
                case 4
                    PercVec = GRD.NodeSelect.percvec; 
                    nPerc = numel(PercVec);
                    TransNodeMeanPerf = zeros(nPerc,1);
                    TransNodeSDPerf  = zeros(nPerc,1);
            end
        end
       
        for PercStep = 1:nPerc

            BinaryGridSelection = ...
                nk_GridOptimumSelector(GD.TR(:,curclass,curlabel), GD.TS(:,curclass,curlabel), ...
                Regul(:,curclass), MD(:,curlabel), GD.FEAT(:,curlabel), GD.Weights(:,curlabel), ...
                GD.DS(:,curlabel), GD.DV(:,curlabel), Ps{curclass}, Pdesc{curclass}, combcell, act, PercVec(PercStep));
            
            if nPerc > 1
               % Concatenate prediction array 
                TransNodeArray = nk_CatNodes(BinaryGridSelection.bestCV1TsPred,curclass);
                [TransNodeMeanPerf(PercStep), TransNodeSDPerf(PercStep)] = nk_ComputeTransNodePerformance(TransNodeArray, label, f, d, curclass);
                fprintf('\nTransnode CV1 ensemble performance at %g%%-percentile: %1.2f (%1.2f)',PercVec(PercStep), TransNodeMeanPerf(PercStep), TransNodeSDPerf(PercStep))
                if feval(op,TransNodeMeanPerf(PercStep),PerfCritBin)
                    GD.BinaryGridSelection{curclass}{curlabel} = BinaryGridSelection;
                    PerfCritBin = TransNodeMeanPerf(PercStep);
                    OptPerc = PercVec(PercStep);
                end
            else
                GD.BinaryGridSelection{curclass}{curlabel} = BinaryGridSelection;
                if isfield(GRD,'NodeSelect')  && GRD.NodeSelect.mode == 3
                    TransNodeArray = nk_CatNodes(BinaryGridSelection.bestCV1TsPred,curclass);
                    [ix,jx] = size(TransNodeArray);
                    % Loop through TransNodeArray (CV1 partitions)
                    for ti = 1 : ix
                        for tj = 1 : jx

                        end
                    end
                end
            end

        end

        if nPerc > 1
           GD.BinaryGridSelection{curclass}{curlabel}.TransNodeMeanPerf = TransNodeMeanPerf;
           GD.BinaryGridSelection{curclass}{curlabel}.TransNodeSDPerf = TransNodeSDPerf;
           switch MODEFL
               case 'classification'
                   if RAND.Decompose ~= 9
                       fprintf('\nClassifier #%g: ',curclass);
                   else
                       fprintf('\nMulti-group classifier:')
                   end
               otherwise
                   fprintf('\nRegressor: ')
           end
           cprintf('black*','Final transnode CV1 ensemble performance %g at %g-percentile using %g node(s).', ...
               PerfCritBin, OptPerc, GD.BinaryGridSelection{curclass}{curlabel}.Nodes);
        end

        if  ndims(GD.BinaryGridSelection{curclass}{curlabel}.bestfeats{1}) > 2
            for nodes=1:GD.BinaryGridSelection{curclass}{curlabel}.Nodes
                GD.BinaryGridSelection{curclass}{curlabel}.bestfeats{nodes} = ...
                    GD.BinaryGridSelection{curclass}{curlabel}.bestfeats{nodes}(:,:,curclass);
            end
        end
        if ndims(GD.BinaryGridSelection{curclass}{curlabel}.bestweights{1}) > 2
            for nodes=1:GD.BinaryGridSelection{curclass}{curlabel}.Nodes
                if ~isempty(GD.BinaryGridSelection{curclass}{curlabel}.bestweights)
                    GD.BinaryGridSelection{curclass}{curlabel}.bestweights{nodes} = ...
                        GD.BinaryGridSelection{curclass}{curlabel}.bestweights{nodes}(:,:,curclass);
                end
            end
        end
        if ndims(GD.BinaryGridSelection{curclass}{curlabel}.bestpred{1}) > 2
            for nodes=1:GD.BinaryGridSelection{curclass}{curlabel}.Nodes
                GD.BinaryGridSelection{curclass}{curlabel}.bestpred{nodes} = ...
                    GD.BinaryGridSelection{curclass}{curlabel}.bestpred{nodes}(:,:,curclass);
            end
        end

        if SAV.savemodel && ndims(GD.BinaryGridSelection{curclass}{curlabel}.bestmodel{1}) > 2
            for nodes=1:GD.BinaryGridSelection{curclass}{curlabel}.Nodes
                GD.BinaryGridSelection{curclass}{curlabel}.bestmodel{nodes} = ...
                    GD.BinaryGridSelection{curclass}{curlabel}.bestmodel{nodes}(:,:,curclass);
            end
        end  

    end

    %% MULTI-GROUP PREDICTORS:
    if MULTI.flag
       fprintf('\n\n**** Grid selection process for multi-group classifier ****')
       if MULTI.BinBind
           fprintf('\nGenerate multi-group predictor from binary classifiers'' optima.')
           [GD.MultiGroupGridSelection{curlabel}, MultiBinBind] = nk_MultiBinBind(GD, label(:,curlabel), f, d, curlabel);
       else
           fprintf('\nUse multi-group grid optimum for parameter selection.')
           % Optionally, apply penalization term on the grid search to avoid
            % overfitting and promote sparse and diverse models

           if isfield(GRD,'OptRegul') && GRD.OptRegul.flag && isfield(GRD.OptRegul,'type')

               switch GRD.OptRegul.type
                   case 1
                       %Regul = mean(GD.C,3);
                       Regul = nk_ComputeRegulFunction(mean(GD.C(:,:,curlabel),2), 1, GRD.OptRegul.RegulTypeComplexity); 
                   case 2
                       % Model diversity
                       if RFE.Filter.SubSpaceFlag && RFE.Filter.SubSpaceStrategy > 1
                           Regul = nk_ComputeRegulFunction(GD.MultiM_DivV(:,curlabel), 1, GRD.OptRegul.RegulTypeDiversity);            
                       else
                           Regul = nk_ComputeRegulFunction(GD.CV2Div(:,:,curlabel), 1, GRD.OptRegul.RegulTypeDiversity); 
                       end
                   case 3
                       RegulC = nk_ComputeRegulFunction(mean(GD.C(:,:,curlabel),2), 1, GRD.OptRegul.RegulTypeComplexity); 
                       if RFE.Filter.SubSpaceFlag > 1 && RFE.Filter.SubSpaceStrategy > 1
                            RegulDiv = nk_ComputeRegulFunction(GD.MultiM_DivV(:,curlabel), 1, GRD.OptRegul.RegulTypeDiversity); 
                       else
                            RegulDiv = nk_ComputeRegulFunction(GD.CV2Div(:,:,curlabel), 1, GRD.OptRegul.RegulTypeDiversity); 
                       end
                       Regul = (RegulDiv + RegulC) / 2;
               end
           else
               Regul = mean(GD.C,2);
           end
           % Now select optimal parameters for the multi-group learning machine
           if nPerc > 1 ; 
               TransNodeMeanPerf = zeros(nPerc,1);
               TransNodeSDPerf  = zeros(nPerc,1);
           end

           for PercStep = 1:nPerc

               MultiGroupGridSelection = ...
                    nk_GridOptimumSelector(GD.MultiTR(:,curlabel), GD.MultiTS(:,curlabel), ...
                            Regul, MD(:,curlabel), GD.FEAT(:,curlabel), GD.Weights(:,curlabel), ...
                            GD.MultiCV2Pred(:,curlabel), GD.MultiCV1CVPred(:,curlabel), ...
                            Ps, Pdesc, combcell, act, PercVec(PercStep));
                
               MultiGroupGridSelection.bestprob    = GD.MultiCV2Prob(MultiGroupGridSelection.Npos,curlabel);         
               MultiGroupGridSelection.bestCV2pred = GD.MultiCV1Pred(MultiGroupGridSelection.Npos,curlabel);
               MultiGroupGridSelection.bestCV2prob = GD.MultiCV1Prob(MultiGroupGridSelection.Npos,:,curlabel);
    
               if nPerc > 1
                   % Concatenate prediction array 
                    TransNodeArray = nk_CatNodes(GD.MultiCV1CVPred(MultiGroupGridSelection.Npos));
                    [TransNodeMeanPerf(PercStep), TransNodeSDPerf(PercStep)] = nk_ComputeTransNodePerformance(TransNodeArray, label, f, d);
                    if feval(op,TransNodeMeanPerf(PercStep),PerfCritMult)
                        GD.MultiGroupGridSelection{curlabel} = MultiGroupGridSelection;
                        PerfCritMult = TransNodeMeanPerf(PercStep);
                        OptPerc = PercVec(PercStep);
                    end
                else
                    GD.MultiGroupGridSelection{curlabel} = MultiGroupGridSelection;
                end

           end
       end

       if nPerc > 1
          GD.MultiGroupGridSelection{curlabel}.TransNodeMeanPerf = TransNodeMeanPerf;
          GD.MultiGroupGridSelection{curlabel}.TransNodeSDPerf = TransNodeSDPerf;
          fprintf('\nFinal CV1-performance %g at %g-percentile using %g node(s).', ...
           PerfCritMult, OptPerc, GD.BinaryGridSelection{curclass}.Nodes)
       end
   
    end
end
          