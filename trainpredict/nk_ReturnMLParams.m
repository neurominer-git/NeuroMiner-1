function [Params, Params_desc] = nk_ReturnMLParams(SVM, GRD, nclass, strout)
global VERBOSE

PX = nk_ReturnParamChain(GRD);
Params = cell(nclass,1);
Params_desc = cell(nclass,1);

for i=1:nclass

    switch SVM.prog 

        case {'MikRVM','MKLRVM','GLMFIT','LIKNON','kNNMEX', 'BLOREG', 'LSTSVM'}

            if VERBOSE, fprintf('\n%s #%g: Slack parameters will be omitted.',strout,i); end

        otherwise

            % Check classifier type (L1/L2 or Nu / One-class)
            ctype = nk_GetLIBSVMClassType(SVM);
            % Check if regression parameters are needed 
            rtype = nk_GetLIBSVMRegrType(SVM);
            switch rtype
                case 1
                    [Params{i}{1}, Params_desc{i}{1} ] = nk_ReturnParam('ML-Eps-SVR parameter(s)',PX.Params_desc, PX.Params);
                case 2
                    [Params{i}{1}, Params_desc{i}{1} ] = nk_ReturnParam('ML-Nu-SVR parameter(s)',PX.Params_desc, PX.Params);
            end
            % This is the slack / nu-SVC parameter of the SVM
            switch ctype 
                case 1
                    [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Nu-SVC parameter(s)',PX.Params_desc, PX.Params);
                otherwise
                    [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Slack/Regularization parameter(s)',PX.Params_desc, PX.Params);
            end
    end
    
    switch SVM.kernel.kernstr
        
        case {' -t 0',' -t 4',' -t 5', 'lin','linear','lin_kernel','none','lin_elm'}

            if VERBOSE, fprintf('\n%s #%g: Kernel parameters will be omitted.',strout,i); end

        case {' -t 1', 'poly', 'polynomial', 'Polynomial', 'polyN', 'hpolyN'}
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Kernel parameter(s)',PX.Params_desc, PX.Params);
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Polynomial degree',PX.Params_desc, PX.Params);
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Polynomial coefficients',PX.Params_desc, PX.Params);
            
        case ' -t 3'
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Kernel parameter(s)',PX.Params_desc, PX.Params);
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Sigmoid coefficients',PX.Params_desc, PX.Params);

        otherwise
            % This is the gamma exponent parameter of the RBF kernel
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Kernel parameter(s)',PX.Params_desc, PX.Params);
    end

    switch SVM.prog
        case 'MEXELM'
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-No. hidden neurons',PX.Params_desc, PX.Params);
        case 'kNNMEX'
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Number(s) of nearest neighbors',PX.Params_desc, PX.Params);
        case 'BLOREG'
            [Params{i}{end+1}, Params_desc{i}{end+1} ] = nk_ReturnParam('ML-Tolerance parameter(s)',PX.Params_desc, PX.Params);
    end

end