function cmdstr = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmdstr)
global LIBSVMTRAIN

if isfield(SVM.(SVM.prog),'Weighting') && SVM.(SVM.prog).Weighting   
    switch MODEFL
        case 'classification'
            % Works only for C-SVC (not nu-SVC) as +C / -C slacks are treated separately
            npos = sum(label == 1); nneg = numel(label) - npos; bal = npos / nneg;
            cmdstr = sprintf('%s -w-1 %g', cmdstr, bal^CMDSTR.WeightFact);             
    end
end
if ~strcmp(LIBSVMTRAIN,'svmtrain291')
    cmdstr = [ cmdstr ' -q 1' ];
end
