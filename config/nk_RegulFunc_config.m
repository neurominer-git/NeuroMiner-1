function RegulFunc = nk_RegulFunc_config(defaultfl)

if ~exist('defaultfl','var') || isempty(defaultfl) 
    defaultfl = false;
end

if ~defaultfl
    RegulType = nk_input('Regularization function',0,'m', ...
        ['none|', ...
         '1 - values|', ...
         'cosinus|', ...
         'sinus|', ...
         '1 - cosinus|', ...
         '1 - sinus|', ...
         'sinus * cosinus|', ...
         '1 - (sinus * cosinus)|', ...
         '(sinus + cosinus) / 2|', ...
         '1 - (sinus + cosinus) / 2'], 1:10);

    switch RegulType
        case 1
            RegulFunc = 'none';
        case 2
            RegulFunc = 'invert';
        case 3
            RegulFunc = 'cos';
        case 4
            RegulFunc = 'sin';
        case 5
            RegulFunc = 'cos_invert';
        case 6
            RegulFunc = 'sin_invert';
        case 7
            RegulFunc = 'sin*cos';
        case 8
            RegulFunc = 'sin*cos_invert';
        case 9
            RegulFunc = 'sin+cos/2';
        case 10
            RegulFunc = 'sin+cos/2_invert';
    end
else
    RegulFunc = 'none';
end

end