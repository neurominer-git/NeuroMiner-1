function P = nk_CreateSVMParamArray(Param, P)

P.SVM                           = Param;
ctype                           = nk_GetLIBSVMClassType(P.SVM);
switch ctype
    case {0,4}
        P.Params{1}             = Param.SlackParam;
        P.Params_desc           = {'SlackParam'};
    case 1
        P.Params{1}             = Param.SlackParam;
        P.Params_desc           = {'NuCParam'};
end
rtype                           = nk_GetLIBSVMRegrType(P.SVM);
if rtype > 0
    if isfield(Param,'NuParam') && strcmp(P.SVM.prog,'LIBSVM')
        P.Params{end+1}         = Param.NuParam;
        P.Params_desc{end+1}    = 'NuParam';
    elseif isfield(Param,'EpsParam')
        P.Params{end+1}         = Param.EpsParam;
        P.Params_desc{end+1}    = 'EpsParam';
    end
end
P.Comb = allcomb(P.Params,'matlab');

end