function SVM = GetSVMStruct(FEATSEL, nclass, curclass)

if iscell(FEATSEL.SVM)
    if numel(FEATSEL.SVM)<nclass
        SVM = FEATSEL.SVM{end};
    else
        SVM = FEATSEL.SVM{curclass};
    end
else
    SVM = FEATSEL.SVM;
end


end