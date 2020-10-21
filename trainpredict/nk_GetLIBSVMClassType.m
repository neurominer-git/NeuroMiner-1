function ctype = nk_GetLIBSVMClassType(SVM)

ctype = 0;

if strcmp(SVM.prog,'LIBSVM')

    switch SVM.LIBSVM.LIBSVMver
        
        case 3            
            switch SVM.LIBSVM.classifier
                case {0,1} %L1 or L2-regularized SVC
                    ctype = 0;
                case 2 % Nu-SVC
                    ctype = 1;
                otherwise % One-class SVC
                    ctype = 2;  
            end
                    
        otherwise
            ctype = SVM.LIBSVM.classifier;

    end
    
elseif strcmp(SVM.prog,'LIBLIN')
    
    ctype = 0;
        
end