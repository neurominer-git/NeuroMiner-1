function [OOCVres, cases] = nk_ExtractPredictions(dat, analysis, dim)

flgcv   = nk_input('Display the CV data, too?',0,'yes|no',[1,0],1);
labels  = dat.label_oocv;
cases   = dat.cases_oocv;
if flgcv
    
    if isfield(analysis{dim}.OOCV,'BinResults')
        frmwrk              = 1;
        nclass              = length(analysis.params.cv.class{1,1});
        R                   = analysis.OOCV.BinResults;
        [ CV2Pred, CV2Ind ] = GetPredictions(analysis, dim, nclass);

    elseif isfield(analysis{dim}.OOCV,'MultiResults')
        frmwrk      = 2;
        nclass      = length(analysis.params.cv.class{1,1});
        R           = analysis.OOCV.MultiResults;
        [ CV2Pred, CV2Ind, MultiCV2Pred ] = GetPredictions(analysis, dim, nclass);
        
    else
        frmwrk      = 3;
        nclass      = 1;
        R           = analysis.OOCV.RegrResults;
        CV2Pred{1}  = analysis.GDdims{dim}.Regr.mean_predictions;
    end

    flgsub = nk_input('Only a specific subgroup in the CV data?',0,'yes|no',[1,0],1);
    if flgsub,
        switch frmwrk
            case {1,2}
                ind = CV2Ind{indclass}(subind);
            case 3
                subind = nk_input('Specifiy index to subgroup observations',0,'e');
                ind = analysis.GDdims{dim}.Regr.index_predictions(subind);
        end
    else
        ind = analysis.GDdims{dim}.Regr.index_predictions;
    end
    fname = 'OOCV & CV prediction results';
    cases = [ cases ; dat.cases(ind) ] ;
    labels = [ labels ; dat.label(ind) ] ;
    R.MeanCV2PredictedValues = [ R.MeanCV2PredictedValues ; analysis.GDdims{dim}.Regr.mean_predictions(ind) ] ;
    flgadd = nk_input('Add the CV2 data to some of the OOCV groups',0,'yes|no',[1,0]);
    if flgadd
        fprintf('\n\nAvailable OOCV groups:')
        fprintf('\n======================')
        for u = 1 : numel(R.Group)
            fprintf('\n%g) %s',u,char(R.Group{u}.GroupName));
        end
        addind = nk_input('Select group index',0,'e');
        R.Group{addind}.ObservedValues = [R.Group{addind}.ObservedValues; dat.label(ind)]; 
        R.Group{addind}.MeanCV2PredictedValues = [R.Group{addind}.MeanCV2PredictedValues; analysis.GDdims{dim}.Regr.mean_predictions(ind)];
        R.Group{addind}.CorrPredictObserved = CC(R.Group{addind}.ObservedValues, R.Group{addind}.MeanCV2PredictedValues);
    else
        R.Group{end+1}.ObservedValues = dat.label(ind); 
        R.Group{end}.MeanCV2PredictedValues = analysis.GDdims{dim}.Regr.mean_predictions(ind);
        R.Group{end}.CorrPredictObserved = CC(R.Group{end}.ObservedValues, R.Group{end}.MeanCV2PredictedValues);
        R.Group{end}.GroupName = dat.groupnames;
    end
else
    fname = 'OOCV prediction results';
end

end

function [ CV2Pred, CV2Ind, MultiCV2Pred ] = GetPredictions(analysis, dim, nclass)

MultiCV2Pred    = [];
CV2Pred         = cell(nclass,1);
CV2Ind          = cell(nclass,1);

for curclass = 1:nclass
    CV2Pred{curclass} = analysis.GDdims{dim}.BinClass{curclass}.mean_predictions;
    CV2Ind{curclass} = analysis.GDdims{dim}.BinClass{curclass}.index_predictions;
end

if isfield(analysis.GDdims{dim},'MultiClass')
    MultiCV2Pred = analysis.GDdims{dim}.MultiClass.multi_probabilitiesCV2;
end

end