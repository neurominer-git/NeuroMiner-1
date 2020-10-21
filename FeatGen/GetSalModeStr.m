function [S, Salgo] = GetSalModeStr(FEATSEL)

Salgo = 0;
switch FEATSEL.salthreshmode 
    case 0 
        switch FEATSEL.salthreshalgo 
            case 1
                S = 'IMRelief';
            case 2
                S = 'kNN';
            case 3
                S = 'RVM';
            case 4
                S = 'SVM';
            case 5
                S = 'LIBLIN';
                
        end
        Salgo = FEATSEL.salthreshalgo;
    case 1
        S = 'FixPercThresh';
    case 2
        S = 'AbsValThresh';
    case 3
        S = 'AdaptiveThresh';    
end          


end
