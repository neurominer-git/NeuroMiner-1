function NumParam = nk_CheckNumParam

global SVM 

NumParam = 0;

if iscell(SVM.kernel.kernstr), kerntyp = SVM.kernel.kernstr{1}; else kerntyp = SVM.kernel.kernstr; end

switch SVM.prog
    case {'LIBSVM','SVMLIT','SVMPRF','CUDSVM'}
        if any(strcmp(kerntyp,{' -t 0',' -t 4',' -t 5', 'lin', 'linear', 'lin_kernel'}))
            NumParam = 1;
        else
            NumParam = 2;
        end
    case {'MKLRVM','MikRVM','MVTRVR'}
        if any(strcmp(kerntyp,{' -t 0',' -t 4',' -t 5', 'lin', 'linear', 'lin_kernel'}))
            NumParam = 0;
        else
            NumParam = 1;
        end

end