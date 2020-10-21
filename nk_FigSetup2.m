function vargout = nk_FigSetup2(fl)

if nargin < 1, fl = 'all'; end
% ***************************** FIGURE SETUP ******************************
switch fl
    
    case 'cv1test'
        vargout{1} = check_cv1test;
    case 'bestcv1test'
        vargout{1} = check_bestcv1test;
    case 'cv2test'
        vargout{1} = check_cv2test;
    case 'bestcv2test'
        vargout{1} = check_bestcv2test;
    case 'mcv1test'
        vargout{1} = check_mcv1test;
    case 'bestmcv1test'
        vargout{1} = check_bestmcv1test;
    case 'mcv2test'
        vargout{1} = check_mcv2test;
    case 'bestmcv2test'
        vargout{1} = check_bestmcv2test;
    case 'err'
    case 'merr'
    case 'learnperf'
        vargout{1} = check_learnperf;
    case 'maxgrid'
        vargout{1} = check_maxgrid;
    case 'all'
        vargout{1} = check_cv1test;
        vargout{2} = check_bestcv1test;
        vargout{3} = check_cv2test;
        vargout{4} = check_bestcv2test;
        vargout{5} = check_mcv1test;
        vargout{6} = check_bestmcv1test;
        vargout{7} = check_mcv2test;
        vargout{8} = check_bestmcv2test;
        vargout{9} = check_learnperf;
        vargout{10} = check_maxgrid;
end

end

function h = check_cv1test

    ho = findobj('Name','CV1 performance across SVM parameters');
    if isempty(ho), 
        h = figure('Name','CV1 performance across SVM parameters','NumberTitle','off'); 
    else
        h = ho;
    end
    
end

function h = check_bestcv1test

    ho = findobj('Name','CV1 performance across SVM parameters (best dim)');
    if isempty(ho), 
        h = figure('Name','CV1 performance across SVM parameters (best dim)','NumberTitle','off'); 
    else
        h = ho;
    end

end

function h = check_cv2test

    ho = findobj('Name','CV2 performance across SVM parameters');
    if isempty(ho), 
        h = figure('Name','CV2 performance across SVM parameters','NumberTitle','off'); 
    else
        h = ho;
    end
    
end

function h = check_bestcv2test

    ho = findobj('Name','CV2 performance across SVM parameters (best dim)');
    if isempty(ho), 
        h = figure('Name','CV2 performance across SVM parameters (best dim)','NumberTitle','off'); 
    else
        h = ho;
    end

end

function h = check_mcv1test
global MULTI

if MULTI.flag
    ho = findobj('Name','CV1 multi-class performance across SVM parameters');
    if isempty(ho), 
        h = figure('Name','CV1 multi-class performance across SVM parameters','NumberTitle','off'); 
    else
        h = ho;
    end
else
    h=[];
end

end

function h = check_mcv2test
global MULTI

if MULTI.flag
    ho = findobj('Name','CV2 multi-class performance across SVM parameters');
    if isempty(ho), 
        h = figure('Name','CV2 multi-class performance across SVM parameters','NumberTitle','off'); 
    else
        h = ho;
    end
else
    h=[];
end
end

function h = check_bestmcv1test
global MULTI

if MULTI.flag
    ho = findobj('Name','CV1 multi-class performance across SVM parameters (best dim)');
    if isempty(ho), 
        h = figure('Name','CV1 multi-class performance across SVM parameters (best dim)','NumberTitle','off'); 
    else
        h = ho;
    end
else
    h=[];
end

end

function h = check_bestmcv2test
global MULTI

if MULTI.flag
    ho = findobj('Name','CV2 multi-class performance across SVM parameters (best dim)');
    if isempty(ho), 
        h = figure('Name','CV2 multi-class performance across SVM parameters (best dim)','NumberTitle','off');
    else
        h = ho;
    end
else
    h=[];
end

end

function hx = check_learnperf
global DR

switch DR.RedMode
    case 'none'
        hx = findobj('Name','Learning performance');
        if isempty(hx), 
            hx = figure('Name','Learning performance','NumberTitle','off'); 
        end;
    otherwise
        hx = findobj('Name','Learning performance (over feature dimensions)');
        if isempty(hx), 
            hx = figure('Name','Learning performance (over feature dimensions)','NumberTitle','off'); 
        end;
end

end

%end
function hg = check_maxgrid
global GRD

if GRD.GridMaxType
    hg = findobj('Name','Max distribution across grid');
    if isempty(hg)
        hg = figure('Name','Max distribution across grid','NumberTitle','off');
    end
else
    hg = [];
end
end