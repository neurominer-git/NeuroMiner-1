function modeflag = nk_DefineModeflag(modeflag)

defmode = 1;
if exist('modeflag','var') && ~isempty(modeflag)
    modes = {'classification','regression'};
    defmode = find(strcmp(modes,modeflag));
    if isempty(defmode),defmode = 1; end
end

modeflag = nk_input('Would you like to perform ...',0,'m', ...
                'classification|regression',[1,2],defmode);
switch modeflag 
    case 1
        modeflag = 'classification';        
    case 2
        modeflag = 'regression';
end
