function [ fl, desc ] = GetMLType(SVM)

fl = [];

if exist('SVM','var') && ~isempty(SVM) && isfield(SVM,'prog')
    
    switch SVM.prog
        case {'MikRVM', 'MKLRVM', 'MVTRVR','FAMRVR'}
            fl = 'RVM'; desc = 'Relevance Vector Machines';
        case {'LIBSVM','SVMLIT','SVMPRF','PROSVM','CUDSVM'}
            fl = 'SVM'; 
            if strcmp(SVM.prog,'LIBSVM')
                switch SVM.LIBSVM.classifier
                    case 0
                        clstr = 'L1-C-SVC'; 
                    case 1
                        switch SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                clstr = 'nu-SVC';
                            case 3
                                clstr = 'L2-C-SVC ';
                        end
                    case 2
                        switch SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                clstr = 'One-Class';
                            case 3
                                clstr = 'nu-SVC';
                        end
                    case 3
                        switch SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                clstr = 'e-SVR';
                            case 3
                                clstr = 'One-Class';
                        end
                    case 4
                        switch SVM.LIBSVM.LIBSVMver
                            case {0,1,2}
                                clstr = 'nu-SVR';
                            case 3
                                clstr = 'e-SVR';
                        end
                    case 5
                        clstr = 'nu-SVR';
                    case 6 
                        clstr = 'L1-SVDD';
                    case 7
                        clstr = 'L2-SVDD';
                end  
                fl = [fl ' (' clstr ')'];
                desc = 'Support Vector Machines (LIBSVM)';
            else
                desc = 'Support Vector Machines';
            end
        case 'LIBLIN'
            switch SVM.LIBLIN.classifier
                case 0
                    clstr = 'L2-reg LR [primal]';
                case 1
                    clstr = 'L2-reg-L2-loss SVC [dual]';
                case 2
                    clstr = 'L2-reg-L2-loss SVC [primal]';
                case 3
                    clstr = 'L2-reg-L1-loss SVC [dual]';
                case 5
                    clstr = 'L1-reg-L2-loss SVC';
                case 6
                    clstr = 'L1-reg LR';
                case 7
                    clstr = 'L2-reg LR [dual]';
                case 11
                    clstr = 'L2-reg-L2-loss eSVR [primal]';
                case 12 
                    clstr = 'L2-reg-L2-loss eSVR [dual]';
                case 13
                    clstr = 'L2-reg-L1-loss eSVR [dual]'; 
            end
            fl = ['LIBLINEAR (' clstr ')']; desc = 'Linear Learners (LIBLINEAR)';
        case 'IMRELF'
            fl = 'IMRelief'; desc ='LOGO algorithm';
        case 'GLMFIT'
            fl = 'GLM'; desc = 'General Linear Model';
        case 'kNNMEX'
            fl = 'kNN'; desc = 'k-Nearest Neighbor';
        case 'MEXELM'
            fl = 'ELM'; desc = 'Extreme Learning Machines';
        case 'BLOREG'
            fl = 'B-LOG-REG'; desc = 'Sparse Bayesian Logistic Regression';
        case 'LSTSVM'
            fl = 'LS-SVM'; desc = 'Least-Squares SVM';
        case 'matLRN'
            fl = sprintf('matLearn [ %s ]',SVM.matLRN.algo{1}); desc ='matLearn library algorithms';
        case 'DECTRE'
            fl = 'Decision Tree'; desc = fl;
        case 'RNDFOR'
            fl = 'Random Forest'; desc = fl;
        case 'GLMNET'
            fl = 'Elastic net'; desc = fl;
        case 'GRDBST'
            fl = 'Gradient boosting'; desc = fl;
        case 'SEQOPT'
            fl = 'Sequence optimizer'; desc = fl;
        otherwise
            fl = SVM.prog; desc ='other algorithm';
    end

end