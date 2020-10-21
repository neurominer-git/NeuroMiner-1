function [ Tr, CV, Ts, Ocv, TrainedParam ] = nk_ReturnAtOptPos(oTr, oCV, oTs, oOcv, Pnt, z)

Ocv = [];
if ~isempty(Pnt.data_ind)
    Ix = Pnt.data_ind(z);
else
    Ix = 1;
end
if iscell(oTr)
    Tr = oTr{Ix}; CV = oCV{Ix}; Ts = oTs{Ix};   
    if ~isempty(oOcv), Ocv = oOcv{Ix}; end
else
    Tr = oTr; CV = oCV; Ts = oTs;   
    if ~isempty(oOcv), Ocv = oOcv; end
end

if nargout == 5
    if ischar(Pnt.TrainedParam) 
        if exist(Pnt.TrainedParam,'file')
            fprintf('\tLoading parameter file.')
            load( Pnt.TrainedParam );
        else
            error('Parameter file %s could not be found.',Pnt.TrainedParam );
        end
    else
        oTrainedParam = Pnt.TrainedParam;
    end
    if ~isempty(Pnt.nA)
        TrainedParam = cell(1,Pnt.nA);
        for a = 1:Pnt.nA
            if isstruct(Pnt.TrainedParam{a})
                TrainedParam{a} = oTrainedParam{a};
            else
                TrainedParam{a} = oTrainedParam{a}{Pnt.train_ind(z,a)};
            end
        end
    else
        TrainedParam = oTrainedParam;
    end
end