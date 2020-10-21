function passed = nk_DefineCaseListIntegrity(cases, tcases)

passed = 1;

if ~isempty(cases)

    fprintf('\n*** Case list integrity check ***')
    if ~isequal(numel(cases),numel(tcases))
        error('Case numbers in current variate differ from previous variates.')
    end
    if ~isequal(cases,tcases)
        for i=1:numel(tcases)
            if ~isequal(cases{i}, tcases{i}), a = '*'; else a = ''; end
            fprintf('\n'); cprintf('red','old: %s, new: %s\t%s ', cases{i}, tcases{i},a)
        end
        fprintf('\n=====================================================================\n');
        warning('Existing case identifiers are not identical to current ones.')
        fprintf('\n\nYou should know EXACTLY what you do right now !!!')
        abortflag = spm_input('Actions',0,'Ignore and Continue|Abort',[0,1]);
        if abortflag, passed = 0; return; end
    else
        cprintf('green','\npassed.')
    end
end