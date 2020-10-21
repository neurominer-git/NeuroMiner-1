function cParams = nk_PrepMLParams(Params, Params_desc, i)
global SVM GRD CMDSTR

% Pull parameters according to the algorithm selected
switch SVM.prog
    case 'matLRN'
        cParams.val = Params(i,:); 
        cParams.desc = Params_desc;
    case 'GLMNET'
        % Use only the first 5 params
        cParams.val = Params(i,1:5); 
        cParams.desc = Params_desc(1:5);
    case 'GRDBST'
        cParams.val = Params(i,1:5); 
        cParams.desc = Params_desc(1:5);
    case 'ROBSVM'
        cParams.val = Params(i,:); 
        cParams.desc = Params_desc;
        CMDSTR.cmd = nk_GenRobSVMCmd(GRD.ROBSVM,cParams);
        CMDSTR.quiet = ' -q'; 
    case 'SEQOPT'
        cParams.val = Params(i,:);
        cParams.desc = Params_desc;
    case 'WBLCOX'
        cParams.val = Params(i,:);
        cParams.desc = Params_desc;
        
   otherwise
        cParams = Params(i,:); 
        switch SVM.prog
            case {'MKLRVM'}
                 %if nvar >1, cParams = repmat(cParams(curclass),1,nvar);end
            case 'MikRVM'
                cParams = nk_ReturnParam('Kernel', Params_desc, Params(i,:));
            case {'LIBSVM','LIBLIN','CCSSVM'}
                 % Convert parameters to char array
                 cParams = num2str(cParams','%1.10f');
                 % Concatenate parameter string
                 cParams = nk_ConcatLIBSVMParamStr(cParams);
                 if SVM.(SVM.prog).Weighting
                    CMDSTR.WeightFact = nk_ReturnParam('Weight Factor', Params_desc, Params(i,:));
                 end
                 if isfield(SVM.(SVM.prog),'MakeInsensitive') && SVM.(SVM.prog).MakeInsensitive
                     CMDSTR.CCLambda = nk_ReturnParam('CC-Lambda', Params_desc, Params(i,:));
                 end
        end
end
