function linsvmfl = determine_linsvm_flag(SVM)

linsvmfl = any(strcmp({'LIBSVM','LIBLIN'}, SVM.prog)) && any(strcmp({' -t 0',' -t 4',' -t 5','lin', 'linear'}, SVM.kernel.kernstr));

end