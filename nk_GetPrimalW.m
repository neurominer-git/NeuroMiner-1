function [w, b] = nk_GetPrimalW(model)
global SVM MODEFL

switch SVM.prog
    
    case 'LIBSVM'
        
        w = model.SVs' * model.sv_coef;
        b = -model.rho;

        if strcmp(MODEFL,'classification')
            if model.Label(1) == -1
              w = -w;
              b = -b;
            end
        end
        w = w';
        
    case 'MikRVM'
        
    case 'LIBLIN'
        w = model.w;
    case 'MSTOOL'
        w = model.w(2:end)';

end