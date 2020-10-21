function [ D, datatype, brainmask, badcoords, threshval, threshop] = getD(FUSEFLAG, inp, n)

if isfield(inp,'stacking') && inp.stacking
    datatype    = 0;
    brainmask   = [];
    badcoords   = false(1,inp.nD);
    threshval   = []; 
    threshop    = [];
    D           = inp.nD;
else
    switch FUSEFLAG
        case {0,1,3}

            switch FUSEFLAG
                case {0, 3}
                    D = inp.X(1).dimvecx(end);
                    datatype    = inp.X.datatype;
                    brainmask   = inp.X.brainmask{1};
                    badcoords   = inp.X.badcoords{1};
                    threshval   = inp.X.threshval; 
                    threshop    = inp.X.threshop;
                case 1
                    D = inp.X(1).dimsizes(n);
                    datatype    = inp.X.datatype(n);
                    brainmask   = inp.X.brainmask{n};
                    badcoords   = inp.X.badcoords{n};
                    threshval   = inp.X.threshval{n}; 
                    threshop    = inp.X.threshop{n};
            end
        case 2
            D = inp.X(n).dimvecx(end);
            datatype    = inp.X(n).datatype;
            brainmask   = inp.X(n).brainmask{1};
            badcoords   = inp.X(n).badcoords{1};
            threshval   = inp.X(n).threshval; 
            threshop    = inp.X(n).threshop;
    end
end