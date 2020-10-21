function cmdstr = nk_ConcatLIBSVMParamStr(Params)

global CMDSTR

%% Build parameter string
cmdstr = CMDSTR.simplemodel;
for i = 1:numel(CMDSTR.ParamStr)
    cmdstr = [ ' -' CMDSTR.ParamStr{i} ' ' strtrim(Params(i,:)) cmdstr ];
end


end