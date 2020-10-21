function [ suggest best_acc best_sv ] = nk_SA (Y, inp, epsilon, init_T, end_T, dt, lambda, big_gamma, also_plot )
% 
%   function [ suggest best_acc best_sv ] =  nk_SA (Y, inp, epsilon, init_T, end_T, dt, lambda, big_gamma, also_plot )
%       (n-fold cross-validation version, with simulated annealing rather than grid search)
% 
% 
% Copyright (C) 2006, Matthew D. Boardman and Thomas P. Trappenberg, Dalhousie University, Nova Scotia, Canada
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
% 

 
global SVM RFE SAV GRD MULTI 

% be well-behaved in a crash:  set output variables to some value
        suggest=0;
        best_acc=0;
        best_sv=0;

    if nargin < 3

        fprintf('\nUsage:\n\n   [ suggest best_acc best_sv ] = sa_svc ( svm, w, c, lcr, lgr, epsilon, n_fold, init_T, end_t, dt, lambda, big_gamma, also_plot )\n\n');

        fprintf('    Optimize cost and gamma parameters for support vector classifier with RBF kernel.\n');
        fprintf('    Uses simulated annealing, extrinsic regularization and intensity-weighted centre of mass.\n');
        fprintf('    This version uses N-fold cross-validation.\n\n');
        
        fprintf('These parameters must be specified:\n');
        fprintf('        svm  = 1 for LIBSVM or 2 for SVM-Light\n');
        fprintf('        w    = a matrix of training data (rows = training samples, columns = input dimensions)\n');
        fprintf('        c    = a vector of training labels (rows = training samples)\n');
        fprintf('\nThese parameters are optional:\n');
        fprintf('        lcr  = log cost range  (default is [ -5.5:15.5] )\n');
        fprintf('        lgr  = log gamma range (default is [-15.5: 5.5] )\n');
        fprintf('    epsilon  = SVM stopping criteria (default is 0.001)\n');
        fprintf('     n_fold  = number of data partions in N-fold cross-validation (default is 5)\n');
        fprintf('     init_T  = initial temperature for simulated annealing (default is 100)\n');
        fprintf('      end_T  = ending temperature for simulated annealing (default is 1)\n');
        fprintf('        dt   = cooling rate (default is 0.1:  must be between 0 and 1!)\n');
        fprintf('     lambda  = weight of support vector ratio in score (default is 1 = equal weight)\n');
        fprintf('   big_gamma = non-linearity for support vector ratio (default is .5 = sqrt, 1 is linear, 2 is square etc.)\n');
        fprintf('   also_plot = 1 to plot the simulated annealing path as you go (default is no)\n\n');

        fprintf('The parameter you will most want to play with is dt, which controls the length of the cooling cycle.\n');
        fprintf('A value of dt=0.1 gives 440 evaluations, whereas dt=0.01 gives 4590 evaluations.  More is better.\n');
        
        fprintf('Copyright (C) 2006, Matthew D. Boardman and Thomas P. Trappenberg, Dalhousie University, Nova Scotia, Canada.\n');
        fprintf('These MATLAB scripts come with ABSOLUTELY NO WARRANTY; for details, see GNU General Public License (license.txt).\n');
        fprintf('This is free software, and you are welcome to redistribute it under certain conditions; see license.txt for details.\n');
        fprintf('\n\n\n');
        
        return;
    end;
    
    ujumps = 0.1;        % accept (small) upward jumps 10% of the time
    reset_jumps = 0.01;  % occasional reset to best positions
    

k = inp.k; kexp = inp.kexp;
g = inp.g; gexp = inp.gexp;
nclass          = inp.nclass;

% Free memory, otherwise program will crash with huge data sets

% Initialize variables
[m,n] = size(k{1});
if n > 1, 
    % in this case different grids will be used for the CV2
    % partitions, given a greater flexibility to the model construction
    % process
    kl = size(k{1},2);
    gl = size(g{1},2);
else
    % in this case one grid will be used for all CV2 partitions, giving
    % greater parameter stability to the parameter selection process
    kl = size(k{1},1);
    gl = size(g{1},1);
end

npos = kl * gl;

% Setup CV2 container variables:
[ix, jx, nclassY] = size(Y);

if nclass ~= nclassY
    error('Number of binary comparisons is not correct. Check either Y or inp');
end

if nargin<5
    GridAct=nk_CVGridSelector(ix,jx);
end

ol = sum(GridAct(:)); % No. of CV2 partitions to evaluate

% ********************* SETUP RESULTS STRUCTURE ***************************
GDanalysis.params.SVM = SVM;
GDanalysis.params.RFE = RFE;
GDanalysis.GridAct = GridAct;
GDanalysis.GDpaths = cell(ol,1);
GDanalysis.C = k;
GDanalysis.C_exp = kexp;
GDanalysis.Gamma = g;
GDanalysis.Gamma_exp = gexp;
GDanalysis.SvmBinComp = nclass;
GDanalysis.NumCV2Part = ol;

GDanalysis.grid.mean_CVPerf = zeros(kl,gl,nclass);
GDanalysis.grid.mean_TSPerf = zeros(kl,gl,nclass);
GDanalysis.grid.mean_TSPerf_JointHypoTargs = zeros(kl,gl,nclass);
GDanalysis.grid.mean_TSPerf_JointHypoDecs = zeros(kl,gl,nclass);
GDanalysis.grid.mean_Err_CVTSPerf = zeros(kl,gl,nclass);

if MULTI.flag
    GDanalysis.grid.MultiCVPerf = zeros(kl,gl);
    GDanalysis.grid.MultiTSPerf = zeros(kl,gl);
    GDanalysis.grid.sdMultiCVPerf = zeros(kl,gl);
    GDanalysis.grid.sdMultiTSPerf = zeros(kl,gl);
end

GDanalysis.bestc = zeros(ix,jx,nclass);
if ~inp.nogexp,GDanalysis.bestg = zeros(ix,jx,nclass); end
GDanalysis.bestTR = cell(nclass,1);
GDanalysis.bestTS = cell(nclass,1);
GDanalysis.bestc = cell(nclass,1);
GDanalysis.bestg = cell(nclass,1);
GDanalysis.bestcexp = cell(nclass,1);
GDanalysis.bestgexp = cell(nclass,1);
for h=1:nclass
    GDanalysis.bestTR{h} = zeros(ix,jx);
    GDanalysis.bestTS{h} = zeros(ix,jx);
    GDanalysis.bestc{h} = zeros(ix*jx,1);
    GDanalysis.bestg{h} = zeros(ix*jx,1);
    GDanalysis.bestcexp{h} = zeros(ix*jx,1);
    GDanalysis.bestgexp{h} = zeros(ix*jx,1);
end

ll = 1;   

for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds
        
        if ~GridAct(f,d), 
            ll = ll +1;
            continue, 
        end;
        
        GD = cell(nclass,1);
        MD = cell(nclass,1);
        
        for h=1:nclass % Loop through binary comparisons

            % ************************ PREPARATIONS ***********************
            fprintf(1,'\n\n*** SVM #%g ***',h);
            
            % CV1 performance for CV2 partition [f,d] on grid [C,Gamma] (or
            % vector [C])
%             GD{h}.TR    = zeros(kl,gl);
%             GD{h}.mMd   = zeros(kl,gl);
%             
%             % CV2 performance for CV2 partition [f,d] on grid [C,Gamma] (or
%             % vector [C])
%             GD{h}.TS    = zeros(kl,gl); 
%             GD{h}.mTS   = zeros(kl,gl); % mean
%             GD{h}.sTS   = zeros(kl,gl); % sd
%             
%             % Generalization error for CV2 partition [f,d] on grid [C,Gamma] (or
%             % vector [C])
%             GD{h}.ERR  = zeros(kl,gl); % mean
            
            % Models
%             MD{h}       = cell(kl,gl);  % models
%             GD{h}.FEAT  = cell(kl,gl);  % selected features for model in MD
%             GD{h}.W2    = cell(kl,gl);  % margin of models in MD
%             GD{h}.Md    = cell(kl,gl);  % Mean distance to hyperplane of models in MD
%             GD{h}.Mm    = cell(kl,gl);  % Normalized margin of models in MD
            
%             % Prediction & Decision values of CV1 and CV2 test data
%             GD{h}.RT    = cell(kl,gl);  % Prediction of CV1 training data
%             GD{h}.DT    = cell(kl,gl);  % Decision values of CV1 training data
%             GD{h}.RV    = cell(kl,gl);  % Prediction of CV1 test data
%             GD{h}.DV    = cell(kl,gl);  % Decision values of CV1 test data
%             GD{h}.RS    = cell(kl,gl);  % Prediction of CV2 test data
%             GD{h}.DS    = cell(kl,gl);  % Decision values of CV2 test data

            C = k{h}; m = size(C,1); if m>1, C = C(ll,:); end
            if inp.nogexp
                G = NaN;
            else
                G = g{h}; m = size(G,1); if m>1, G = G(ll,:); end
            end
            if RFE.Filter.flag
                FilterSubSets = nk_CreateSubSets(Y{f,d,h},RFE.Filter.SubSpaceFlag);
            else
                FilterSubSets = [];
            end
    
            %--------------------------------------------------------------
            % run the cross validation for simulated annealing approach
            %--------------------------------------------------------------

            % set initial point
            T = init_T;
            log_gamma = 0;
            log_cost = 0;
            i=0;  % drop temp every 10 iterations
            j=0;  % number of accuracy points found

            % guess that there will be 3000 accuracy points in the path (more allocated as we go)
            acc = zeros(3000,3);   
            true_acc = zeros(3000,3);
            true_sv = zeros(3000,3);

            current_acc = 0;
            best_acc = -1e6;
            best_true_acc = 0;
            best_true_sv = 0;
            best_num = 1;

            min_lcr = min(C);  % so that we just calculate this once
            max_lcr = max(C);
            min_lgr = min(G);
            max_lgr = max(G);

            g_factor = (max_lgr-min_lgr)/init_T;   % so that we just calculate this once
            c_factor = (max_lcr-min_lcr)/init_T;

            if  also_plot > 0

                figure(1); hold all;
                title('Simulated Annealing Path');
                xlabel('log_e(cost)');
                ylabel('log_e(gamma)');
                axis ([ min_lcr max_lcr min_lgr max_lgr ]); 
            end;


            %-------------------------------MAIN LOOP -------------------------------------------
            while ( T > end_T )

                % determine a new jump point at random
                new_log_gamma = randn()*T*g_factor;
                new_log_cost  = randn()*T*c_factor;

                % introduce a small bias towards the origin (0,0) to ease complexity of final solution
                bias_gamma = (-log_gamma)*T*g_factor/10;
                bias_cost  = (-log_cost)*T*c_factor/10;

                distance = sqrt(new_log_gamma^2 + new_log_cost^2 + bias_gamma^2 + bias_cost^2);
                new_log_gamma = new_log_gamma + log_gamma + bias_gamma;
                new_log_cost = new_log_cost + log_cost + bias_cost;

                % check that we're within the lcr and lgr bounds: else jump randomly
                if new_log_gamma > max_lgr  
                    new_log_gamma = rand()*(max_lgr-min_lgr)+min_lgr;  
                end;
                if  new_log_gamma < min_lgr
                    new_log_gamma = rand()*(max_lgr-min_lgr)+min_lgr;  
                end;
                if new_log_cost > max_lcr 
                    new_log_cost = rand()*(max_lcr-min_lcr)+min_lcr;  
                end;
                if new_log_cost < min_lcr
                    new_log_cost = rand()*(max_lcr-min_lcr)+min_lcr;  
                end;

                % create parameter string for this run

                %if ~strcmp(SVM.prog,'CUDSVM')
                    kstr = num2str(exp(new_log_cost),'%1.10f');
                    gstr = num2str(exp(new_log_gamma),'%1.10f');
                    estr = num2str(epsilon,'%1.10f');

                %end  

        %-----------------------BEGIN TRIAL---------------------------

               % get cross-validation accuracy and number of support vectors
               [CVperf, TSperf] = nk_CVPermFold(Y{f,d,h}, kstr, gstr, FilterSubSets);
                
%                if RFE.Wrapper.flag % Assign data according to feature selection strategy
% 
%                     [GD{h}.TR(j,i), ...
%                         MD{h}{j,i}, ...
%                         GD{h}.FEAT{j,i}, ...
%                         GD{h}.W2{j,i}, ...
%                         GD{h}.Md{j,i}, ...
%                         GD{h}.Mm{j,i}, ...
%                         GD{h}.RT{j,i}, ...
%                         GD{h}.DT{j,i}, ...
%                         GD{h}.RV{j,i}, ...
%                         GD{h}.DV{j,i}, ...
%                         meanvec, stdvec, xaxlb] = ...
%                         gridsearch_helper(CVperf.Wrapper, RFE.Wrapper, RFE.CV2Class.EnsembleStrategy);
% 
%                 elseif RFE.Filter.flag
% 
%                     [GD{h}.TR(j,i), ...
%                         MD{h}{j,i}, ...
%                         GD{h}.FEAT{j,i}, ...
%                         GD{h}.W2{j,i}, ...
%                         GD{h}.Md{j,i}, ...
%                         GD{h}.Mm{j,i}, ...
%                         GD{h}.RT{j,i}, ...
%                         GD{h}.DT{j,i}, ...
%                         GD{h}.RV{j,i}, ...
%                         GD{h}.DV{j,i}, ...
%                         meanvec, stdvec, xaxlb] = ...
%                         gridsearch_helper(CVperf.Filter, RFE.Filter, RFE.CV2Class.EnsembleStrategy);
%                 end

        %         if svm == 1     % LIBSVM
        %             ac = svmtrain(c,w,cat(2,params,s_fold));
        %         else            % SVMLight
        %             ac = svmlearn(w,c,cat(2,params,s_fold));
        %         end;
        % 
        %        % train model (on all training data) to get number of SV
        %         if svm == 1     % LIBSVM
        %             model = svmtrain(c,w,params);
        %         else            % SVMLight
        %             model = svmlearn(w,c,params);
        %         end;

                new_acc = ac;

                new_sv = model.totalSV;        
                %new_sv = model{3};            % if above line causes an error, use this instead (old versions etc.)

                fprintf('Number of support vectors = %d\n', new_sv);

        %----------------------END TRIAL--------------------------------        

                % these may be kept for later display purposes:
                true_acc(j+1,1)=new_log_gamma;
                true_acc(j+1,2)=new_log_cost;
                true_acc(j+1,3)=new_acc;

                true_sv(j+1,1)=new_log_gamma;
                true_sv(j+1,2)=new_log_cost;
                true_sv(j+1,3)=new_sv;

                % calculate score 
                % note: lambda=1 for equal weight
                new_acc = new_acc - lambda*100*((new_sv/size(w,1))^(big_gamma));

                fprintf ('#%4d: (Cost=exp(%3.1f),Gamma=exp(%3.1f),Dist=%.3f): Score=(%.2f->%.2f):Temp=%.1f:', j+1, new_log_cost, new_log_gamma, distance, current_acc, new_acc, T );

                % add to the list of known positions
                j=j+1;
                acc(j,1) = new_log_gamma;      
                acc(j,2) = new_log_cost;
                acc(j,3) = new_acc; 

                % do we want to accept this jump?
                if (new_acc >= current_acc) | (j==1)
                    % descent (or flat):  accept always
                    accept = 1;
                    fprintf('Good jump');
                else
                    % ascent:  maybe accept, if it's not too big of a jump upwards
                    if (abs(current_acc-new_acc) < (current_acc * 0.1) )  & ( rand() < ujumps )  % allow 10% downwards (10% of the time)
                   % if  rand() < ujumps   % allow upward jumps (10% of the time)
                        accept = 1;
                        fprintf('Accepting bad jump!!');
                    else
                        accept = 0;
                        fprintf('Bad jump');
                    end;
                end;

                % plot the jump patterns
                if also_plot > 0

                    jc = 'k';   % jumps not taken
                    figure(1);
                    plot(new_log_cost,new_log_gamma,'o','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1,'MarkerFaceColor','k');

                    if accept > 0 
                        plot([ log_cost new_log_cost ],[ log_gamma new_log_gamma ],jc);
                    end;

                    drawnow;
                end;

                % perform jump if accepted
                if accept > 0
                    current_acc = new_acc;
                    log_gamma = new_log_gamma;
                    log_cost = new_log_cost;    
                end;

               % check best accuracy found so far
                if current_acc > best_acc
                    best_acc = current_acc;
                    clear best_list;
                    best_num = 1;
                    best_list(1,1) = log_cost;
                    best_list(1,2) = log_gamma;
                    best_true_acc = true_acc(j,3);
                    best_true_sv = new_sv;
                elseif current_acc == best_acc
                    best_num = best_num + 1;
                    best_list ( best_num, 1 ) = log_cost;
                    best_list ( best_num, 2 ) = log_gamma;
                end;

                fprintf('\t\t(Best:%.2f=%.2f+%.2f,lG=%.4f,lC=%.4f)\n', best_acc, best_true_acc, best_true_sv, best_list(1,2), best_list(1,1) );

                % drop the temperature by dt every 10 iterations 
                % (regardless of jump acceptance!)
                i=i+1;
                if i >= 10
                    fprintf('*** Dropping temperature from %.1f ', T);
                    T = T*(1-dt);     
                    i = 0;
                    fprintf('to %.1f\n', T);
                end;

                % maybe (but rarely!) reset to a random spot on the best_list,
                %   IFF current accuracy is less than 90% of best accuracy found
                if ( rand() < reset_jumps ) & ( current_acc < best_acc*.9 )

                    jump_to = ceil( rand() * best_num );  % pick a random member of the best_list

                    current_acc = best_acc;
                    log_gamma = best_list ( jump_to, 2 );
                    log_cost = best_list ( jump_to, 1 );

                    fprintf('*** Randomly jumping back to best_list (%d of %d)\n', jump_to, best_num );
                end;

            end;


        %-----------------------END MAIN LOOP--------------------------



        %--------------------------------------------------------------
        % calculate score-weighted centre of mass
        %--------------------------------------------------------------

            fprintf('End of simulated annealing.\n\nCalculating suggested parameters ... ');

            % get suggested point as minimum distance from (0,0)
            for i=1:best_num
                best_list(i,3) = sqrt(best_list(i,1)^2+best_list(i,2)^2);
            end;
            best_list=sortrows(best_list,3);
            suggest=best_list(1,1:2);

            % get overall best suggested point(s)
            clear best_list;
            best_num=0;
            psuggest=suggest;
            suggest=zeros(2,1);
            denom=0;

            acc=acc(1:j,:);   % cut off excess points if needed (since we guessed at the number of iterations)

            for i=1:j               
                 if (sqrt((acc(i,1)-psuggest(1,2))^2+(acc(i,2)-psuggest(1,1))^2) < 1 )   % params are within log radius of 1?
                    if ((best_acc > 0) & (acc(i,3) >= best_acc*.98)) | (best_acc < 0) & (acc(i,3) >= best_acc*1.02) % score is within 2% of best point?

                        if best_acc > 0
                            suggest(1)=suggest(1)+(acc(i,3))*acc(i,2);  % cost
                            suggest(2)=suggest(2)+(acc(i,3))*acc(i,1);  % gamma
                            denom=denom+(acc(i,3));                     % score
                        else    % reverse if scores are negative (not really ideal ... any ideas?)
                            suggest(1)=suggest(1)+(100+acc(i,3))*acc(i,2);  % cost
                            suggest(2)=suggest(2)+(100+acc(i,3))*acc(i,1);  % gamma
                            denom=denom+(100+acc(i,3));                     % score
                        end;

                        best_num=best_num+1;

                        best_list(best_num,1)=acc(i,2);  % cost
                        best_list(best_num,2)=acc(i,1);  % gamma
                        best_list(best_num,3)=acc(i,3);  % score
                    end;
                 end;
            end;

            suggest=suggest/denom;

            % run SVM again with these parameters
            params = sprintf ( '-s 0 -t 2 -g %.8f -c %.8f -e %.8f', exp(suggest(2)), exp(suggest(1)), epsilon );
            s_fold = sprintf ( ' -v %d', n_fold);

            %---------- BEGIN TRIAL -----------

               % get cross-validation accuracy
                if svm == 1     % LIBSVM
                    ac = svmtrain(c,w,cat(2,params,s_fold));
                else            % SVMLight
                    ac = svmlearn(w,c,cat(2,params,s_fold));
                end;

                % train model (on all training data) to get number of SV
                if svm == 1     % LIBSVM
                    model = svmtrain(c,w,params);
                else            % SVMLight
                    model = svmlearn(w,c,params);
                end;

                new_acc = ac;

                new_sv = model.totalSV;        
                %new_sv = model{3};            % if above line causes an error, use this instead (old versions etc.)

            %---------- END TRIAL -----------

            suggest_acc = new_acc;
            suggest_sv = new_sv;

           % keep overall best points too, for display or comparison
            best_acc = best_true_acc;
            best_sv = best_true_sv;

            %--------------------------------------------------------------
            % plot the results as a 3D scatterplot
            %--------------------------------------------------------------

            %  if also_plot > 0
            %     figure(1);
            %     scatter3(true_acc(:,2),true_acc(:,1),true_acc(:,3),8,true_acc(:,3),'filled'); %view(2);
            %     colorbar; grid off; colormap gray;
            %     set(gca,'Color',[ 0 0 .4 ]);
            %     xlabel('log_e(cost)');
            %     ylabel('log_e(gamma)');
            %     title('Points Analyzed (3D scatter plot): Leave-One-Out Accuracy');
            %     axis ([ min_lcr max_lcr min_lgr max_lgr ]); 
            %  end;


            %--------------------------------------------------------------
            % detect failure, otherwise print results
            %--------------------------------------------------------------

            if best_acc == 0
                fprintf ( '\n\n!!! Resulting accuracy could not be measured (same over the entire range of values for cost and gamma).\n');
                return;
            else
                fprintf ('\n\n--------------------\nDone! Top few points:\n');

                for i=1:best_num
                    fprintf ('   Score = %.1f    Gamma = exp(%3.3f) = %.8f\t  Cost = exp(%3.3f) = %.8f\n', ...
                        best_list(i,3), best_list(i,2), exp(best_list(i,2)), best_list(i,1), exp(best_list(i,1)) );
                end;

                fprintf('\nBest point found: \n');
                fprintf ('   Acc = %.1f  SV = %d/%d  Gamma = exp(%3.3f) = %.8f\t  Cost = exp(%3.3f) = %.8f\n', ...
                    best_acc, best_sv, size(w,1), psuggest(1,2), exp(psuggest(1,2)), psuggest(1,1), exp(psuggest(1,1)) );

                fprintf('\nSuggested optimum point: (intensity-weighted centre of mass of best points)\n');
                fprintf ('   Acc = %.1f  SV = %d/%d  Gamma = exp(%3.3f) = %.8f\t  Cost = exp(%3.3f) = %.8f\n\n\n', ...
                    suggest_acc, suggest_sv, size(w,1), suggest(2), exp(suggest(2)), suggest(1), exp(suggest(1)) );

            end;
        end
    end
end
