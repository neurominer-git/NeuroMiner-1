function rtype = nk_GetLIBSVMRegrType(SVM)

rtype = 0;

if strcmp(SVM.prog,'LIBSVM')

    switch SVM.LIBSVM.LIBSVMver
        
        case 3
            
            switch SVM.LIBSVM.classifier
                
                case 4
                    rtype = 1;
                case 5
                    rtype = 2;
                    
            end
                    
        otherwise
            
            switch SVM.LIBSVM.classifier
                
                case 3
                    rtype = 1;
                    
                case 4
                    rtype = 2;
            end

    end
    
elseif strcmp(SVM.prog,'LIBLIN')
    
    switch SVM.LIBLIN.classifier
        
        case {11,12,13}
            rtype=1;
    end


end