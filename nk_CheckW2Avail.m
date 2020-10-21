function w2flag = nk_CheckW2Avail(SVM, MODEFL)

w2flag = false;

switch SVM.prog
    case 'LIBSVM'
        if SVM.LIBSVM.LIBSVMver == 3 && strcmp(SVM.kernel.kernstr, ' -t 0') && ~strcmp(MODEFL,'regression'), w2flag = true; end 
    case 'LIBLIN'
        %w2flag = true;
end

end