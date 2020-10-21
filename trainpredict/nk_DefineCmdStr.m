function vargout = nk_DefineCmdStr(SVM, MODEFL)

vargout = [];

switch SVM.prog
    
    case 'MikRVM'
        vargout = '';
    
    case {'LIBSVM','CCSSVM'}
        if ~SVM.LIBSVM.LIBSVMver
            vargout.quiet = ' -q';
        else
            vargout.quiet = '';
        end
        % Define main parameter string
		vargout.simplemodel = sprintf(' -s %g%s -m %g -e %g -h %g -b %g', ...
                                        SVM.LIBSVM.classifier, ...
                                        SVM.kernel.kernstr, ...
                                        SVM.LIBSVM.Optimization.m, ...
                                        SVM.LIBSVM.Optimization.e, ...
                                        SVM.LIBSVM.Optimization.h, ...
                                        SVM.LIBSVM.Optimization.b);
        
        % Define parameter strings depending on LIBSVM version and
        % classification or regression setup
        switch MODEFL

            case 'regression'
                rtype = nk_GetLIBSVMRegrType(SVM);
                switch rtype
                    case 1 % Epsilon-SVR
                        vargout.ParamStr = {'p','c'};
                    case 2 % Nu-SVR        
                        vargout.ParamStr = {'n','c'};
                end

            case 'classification'
                ctype = nk_GetLIBSVMClassType(SVM);
                switch ctype
                    case 0 % L1- or L2-reg. SVC
                        vargout.ParamStr = {'c'};  
                    case {1,2} % Nu-SVC or one-class SVM
                        vargout.ParamStr = {'n'};                        
                end
        end
        switch SVM.prog
            case 'LIBSVM'
                switch SVM.kernel.kernstr
                    case ' -t 0'
                        vargout.ParamStr = vargout.ParamStr ;
                    case ' -t 1'
                        vargout.ParamStr = [vargout.ParamStr 'g', 'd', 'r'];
                    case ' -t 2'
                        vargout.ParamStr = [vargout.ParamStr 'g'];
                    case ' -t 3'
                         vargout.ParamStr = [vargout.ParamStr 'g', 'r'];
                    otherwise
                        switch SVM.LIBSVM.LIBSVMver
                            case {0,1,2} 
                                error('Precomputed kernels are not supported in this version of NM');

                            case 3  %LIBSVM 2.89 PLUS
                                switch SVM.kernel.kernstr
                                    case {' -t 4', ' -t 5'}
                                        vargout.ParamStr = [vargout.ParamStr 'r'];
                                    case {' -t 6', ' -t 7'}
                                        vargout.ParamStr = [vargout.ParamStr 'g'];
                                    case ' -t 8'
                                        error('Precomputed kernels are not supported in this version of NM');
                                end
                        end
                end
            case 'CCSSVM'
                vargout.ParamStr = [vargout.ParamStr 't 4'];
        end
            
    case 'KPCSVM'
        
        vargout.simplemodel = [SVM.classifier ...
            SVM.kernel.kernstr ...
            SVM.Optimization.m ...
            SVM.Optimization.e ...
            SVM.Optimization.h ];
        
    case 'LIBLIN'
        
                vargout.simplemodel = [ ' -s ' num2str(SVM.LIBLIN.classifier) ...
                                        ' -e ' num2str(SVM.LIBLIN.tolerance) ];
                
        switch MODEFL
            case 'classification'
                vargout.ParamStr = {'c'}; 
            case 'regression'
                vargout.ParamStr = {'p', 'c'}; 
        end
        
        vargout.quiet = ' -q 1';
        
    case 'SVMLIT'
        
        vargout.simplemodel = [SVM.kernel.typ ...
	      SVM.Verbose ... 
	      SVM.BiasHyperplane ...
	      SVM.Optimization.q ...
	      SVM.Optimization.n ...
	      SVM.Optimization.m ...
	      SVM.Optimization.e ...
	      SVM.Optimization.h ...
	      SVM.RemIncons];

        vargout.crossval = [SVM.kernel.typ ...
	      SVM.Verbose ... 
	      SVM.BiasHyperplane ...
	      SVM.Optimization.q ...
	      SVM.Optimization.n ...
	      SVM.Optimization.m ...
	      SVM.Optimization.e ...
	      SVM.Optimization.h ...
	      SVM.RemIncons];
        
    case 'RNDFOR'
        
        switch MODEFL
            case 'classification'
                vargout.RFtrain = str2func('classRF_train');
                vargout.RFtest = str2func('classRF_predict');
            case 'regression'
                vargout.RFtrain = str2func('regRF_train');
                vargout.RFtest = str2func('regRF_predict');
        end
        
end