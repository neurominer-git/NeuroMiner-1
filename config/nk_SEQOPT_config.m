function [act, SEQOPT ] = nk_SEQOPT_config(SEQOPT, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end

if ~defaultsfl 
    
    if ~isfield(SEQOPT,'C'),            SEQOPT.C = NaN; end
    if ~isfield(SEQOPT,'PerfMode'),     SEQOPT.PerfMode = 1; end
    if ~isfield(SEQOPT,'AnchorType'),   SEQOPT.AnchorType = 1; end
    if ~isfield(SEQOPT,'ReplaceMode'),  SEQOPT.ReplaceMode = 1; end
    if ~isfield(SEQOPT,'Mode'),         SEQOPT.Mode = 1; end
    
    MODEARR = {'The population to be propagated','The entire population'}; MODESTR = MODEARR{SEQOPT.Mode};
    MnuStr = sprintf('Optimization based on ... [ %s ]', MODESTR);                                                    MnuAct = 1;
    
    REPLACEARR = {'Replacement','Mean across predictions'};  REPLACESTR = REPLACEARR{SEQOPT.ReplaceMode};
    MnuStr = sprintf('%s|Define computation mode for combining sequential predictions [ %s ]',MnuStr, REPLACESTR); MnuAct = [MnuAct 2];
    
    ANCHORARR = {'Decision boundary','Median'}; ANCHORSTR = ANCHORARR{SEQOPT.AnchorType};
    MnuStr = sprintf('%s|Define propagation algorithm''s anchor [ %s ]',MnuStr, ANCHORSTR);                    MnuAct = [MnuAct 3];
    
    OPTIMARR = {'Optimization criterion','Mean decision distance between classes'}; OPTIMSTR = OPTIMARR{SEQOPT.PerfMode};
    MnuStr = sprintf('%s|Define optimisation mode [ %s ]', MnuStr, OPTIMSTR);                                  MnuAct = [MnuAct 4];
    
    if isnan(SEQOPT.C)
        SEQSTR = 'Sequence undefined';
    else
        SEQSTR = sprintf('%g prediction sequence(s) defined, maximal sequence length: %g predictors',size(SEQOPT.C,1), size(SEQOPT.C,2));
    end
    MnuStr = sprintf('%s|Sequences to be tested [ %s ]', MnuStr, SEQSTR);                                      MnuAct = [MnuAct 5];
    
    nk_PrintLogo
    act = nk_input('Select sequence optimization parameter',0,'mq', MnuStr, MnuAct, 1);
    switch act
        case 1
            SEQOPT.Mode         = nk_input('Define training population flag for optimization',0,'m','The population to be propagated|The entire population',[1,2],SEQOPT.Mode);
        case 2
            SEQOPT.ReplaceMode  = nk_input('Sequential decision computation',0,'m','Replacement|Mean across predictions',[1,2],SEQOPT.ReplaceMode);
        case 3
            SEQOPT.AnchorType   = nk_input('Anchor stepping to model''s decision boundary or decision score median',0,'m','Decision boundary|Median',[1,2], SEQOPT.AnchorType);
        case 4
            SEQOPT.PerfMode     = nk_input('Optimization mode',0,'m','Optimization criterion|Mean decision distance',[1,2], SEQOPT.PerfMode);
        case 5
            SEQOPT.C            = nk_input('Define sequences matrix to be evaluated',0,'e');
    end
else
    act =0; SEQOPT = struct('C', NaN, 'PerfMode', 1, 'AnchorType', 1, 'ReplaceMode', 1, 'Mode', 1); 
end


