function JH = nk_GridSearchEnsembleAdaBoost(Y, GD)
%
% function JH = nk_GridSearchEnsembleAdaBoost(Y, GD, C, G)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Inputs:
% =======
% Y     = Input Data structure (Y.Tr{x,y} = Training Data, Y.TrL{x,y} = Training Data
%                               Labels) produced by nk_Preprocess
% GD    = CVdatamat structure (computed beforehand by NeuroMiner, contains
%                               decisions of trained SVMs, diversity,
%                               performance, complexity info on CV1
%                               training and test data, as well as CV2 test
%                               data
% C     = Slack parameter vector
% G     = Kernel parameter vector
% C & G defined the optimization grid for SVM training
%
% Algorithm:
% ==========
% 
% 1) Define an optimization grid of diversity versus complexity / accuracy
% 2) Loop through this grid and find ensembles that fit into the current
%       grid cell
% 3) Extract the base learners of these ensembles (SVMs with different
%       slack / kernel params, probably) and put their hypotheses into a
%       single matrix
% 4) Forward this matrix to AdaBoost
% 5) In AdaBoost, find optimum weighting of hypotheses on the CV1 training
%       data (CV1 test data ?)
% 6) Use the weights of the optimized committee to predict CV1 test data
% 7) Determine the optimum grid parameters (diversity, complexity/accuracy) 
%       according to the CV1 test data and use the AdaBoost ensemble of
%       this grid cell to predict the CV2 test data
% 8) Return the results in JH.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris & Nikhil Singh 6/2010, BRAINsvm software
% 
global ADA RFE

nclass = size(GD,3);
[iperms,ifolds] = size(Y.Tr);
JH = cell(nclass,1);
ADA.Div = 0.3:0.1:0.8;
ADA.Comp = 0.2:0.1:0.8;
ADA.Perf = 0.5:0.1:0.8;
ADA.MaxIter = 5000;

lDiv = length(ADA.Div); lComp = length(ADA.Comp); lPerf = length(ADA.Perf);

% Loop through binary classifiers
h1 = figure;
h2 = figure;
[lC, lG] = size(GD.FEAT);
if ~exist('GD.DivT','var'), GD.DivT = cell(lC,lG,nclass); DivFl=1; end
for h=1:nclass
    
    JH{h}.final_hyp = cell(lDiv, lPerf);
    JH{h}.weights = cell(lDiv, lPerf);
    JH{h}.distr = cell(lDiv, lPerf);
    JH{h}.err = zeros(lDiv, lPerf);
    JH{h}.CV2perf = zeros(lDiv, lPerf);
    %TrN = size(GD{h}.DT{1,1}{1,1},1);
    %CVN = size(GD{h}.DV{1,1}{1,1},1);
    %TsN = size(GD{h}.DS{1,1}{1,1},1);  
    %TrDx = zeros(TrN,cols); TsDx = zeros(TsN,cols); CVDx =
    %zeros(CVN,cols);
    TsL = Y.TsL{h}; TsInd = Y.TsInd{h};
    %  ************* Aggregate hypotheses and decision values *************
    %  ********************* across parameter grid ************************
    
    % Loop through diversity/complexity grid
    for o=1:lDiv-1
        
        for p=1:lPerf-1
    
            % Loop through Slack, Kernel-Param grid to find ensembles that meet the
            % current diversity and complexity/accuracy window
            TrDx = []; TsDx = []; CVDx = []; TrDiv=[];
            fprintf('\nSearching for ensembles that fit to grid window [div = %g->%g, perf = %g->%g]\n', ...
                ADA.Div(o), ADA.Div(o+1), ADA.Perf(p), ADA.Perf(p+1));
            
            for i=1:lC

                for j=1:lG

                    %fprintf(['\nCheck whether ensemble: SVM #%g, C/Gamma ' ...
                        %'%1.4f; %1.4f fits grid window [div=%g->%g,compl=%g->%g]:\t'], ...
                        %h, i, j, ADA.Div(o), ADA.Div(o+1), ADA.Div(p), ADA.Div(p+1))
                    fprintf('.')
                    if DivFl, GD.DivT{i,j,h} = zeros(iperms,ifolds); end
                        nens=0;
                    for k=1:iperms

                        for l=1:ifolds
                            
                            CVInd = Y.CVInd{k,l}{h};
                            TrInd = Y.TrInd{k,l}{h};
                            
                            TrL = Y.TrL{k,l}{h};
                            
                            CVL = Y.CVL{k,l}{h};
                            
                            % Check for diversity & complexity
                            if DivFl,
                                ET = sign(GD.DT{i,j}{k,l,h}); ET = ET(TrInd,:);
                                EV = sign(GD.DV{i,j}{k,l,h}); EV = EV(CVInd,:);
                                if RFE.ClassRetrain
                                    L = [TrL; CVL];
                                    E = [ET; EV];
                                else
                                    L = TrL;
                                    E = ET;
                                end
                                GD.DivT{i,j,h}(k,l) = nk_DiversityKappa(E, L);
                            end
                            hEPerf = nk_EnsPerf(E, L); if hEPerf > 1, hEPerf=hEPerf/100;end;
                            %fprintf('Diversity = %g, Complexity = %g', GD{h}.DivT{i,j}(k,l), GD{h}.C(i,j));
                            if ~isnan(GD.DivT{i,j,h}(k,l)) && ...
                                    ((GD.DivT{i,j,h}(k,l) < ADA.Div(o) || GD.DivT{i,j,h}(k,l) >= ADA.Div(o+1)) && ...
                                    hEPerf < ADA.Perf(p) || hEPerf >= ADA.Perf(p+1)),
                                continue
                            end
                            fprintf('+')
                            nens=nens+1;
                            
                            trdx = GD.DT{i,j}{k,l,h}(TrInd,:); 
                            cvdx = GD.DV{i,j}{k,l,h}(CVInd,:); 
                            tsdx = GD.DS{i,j}{k,l,h}(TsInd,:);
                            
                            switch RFE.ClassRetrain % 

                                case 0

                                    trl = TrL;
                                    cvl = CVL;
                                    TrDx = [TrDx trdx]; CVDx = [CVDx cvdx]; TsDx = [TsDx tsdx];
                                    
                                case 1
                                    trdx = [trdx; cvdx];
                                    trl = [TrL; CVL];
                                    TrDx = [TrDx trdx]; TsDx = [TsDx tsdx];
                            end
                            %TrDiv = [TrDiv GD{h}.DivT{i,j}(k,l)];
                        end
                    end
                end
            end
            fprintf('\n# of component classifiers in committee: %g', size(TrDx,2))
            
            if size(TrDx,2) > 0
            
                % Invoke AdaBoost
                fprintf('\nApplying AdaBoost on the selected committee.')
                if size(TrDx,2) > ADA.MaxIter
                    TrDx = TrDx(:,1:ADA.MaxIter);
                    TsDx = TsDx(:,1:ADA.MaxIter);
                    if ~RFE.ClassRetrain
                        CVDx = CVDx(:,1:ADA.MaxIter);
                    end
                    %MaxIter = ADA.MaxIter;
                else
                    %MaxIter = size(TrDx,2);
                end
                [JH{h}.final_hyp{o,p}, JH{h}.distr{o,p}, JH{h}.weights{o,p}, JH{h}.err(o,p)] = ns_RealAda(TrDx, trl, [], h1);
                
                % Predict CV2 test data
                % Apply weights to test data hypotheses
                ECV2    =  sign(TsDx).*repmat(JH{h}.weights{o,p},size(TsDx,1),1);
                hECV2   =  sign(sum(ECV2,2));
                JH{h}.CV2perf(o,p) = nk_EnsPerf(hECV2, TsL);
                
            else
                JH{h}.final_hyp{o,p}=[]; JH{h}.distr{o,p}=[]; JH{h}.weights{o,p}=[]; JH{h}.err(o,p) = NaN;
            end
            fprintf('\nCompleted computations for grid cell [div = %g->%g, perf = %g->%g].\n', ...
                ADA.Div(o), ADA.Div(o+1), ADA.Perf(p), ADA.Perf(p+1));
            
            set(0,'CurrentFigure',h2); 
            subplot(2,1,1); contourf(JH{h}.err); 
            colorbar('Location','SouthOutside');
            subplot(2,1,2); contourf(JH{h}.CV2perf);
            colorbar('Location','SouthOutside');
            drawnow 
        end
    end
    
    
end

return
