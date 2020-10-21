function rACTPARAM = nk_ReturnActionParam(ACTPARAM, actstr)

lact = numel(ACTPARAM);
rACTPARAM = [];

for ac=1:lact
    act = ACTPARAM{ac}.cmd;
    if strcmp(act,actstr) 
        rACTPARAM = ACTPARAM{ac};
        break
    end
end

end