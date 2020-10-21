% =========================================================================
% FORMAT [adapthresh, Crit, mx] = nk_AdaptiveThreshold2(Y, labels, V, adap,
% Ynew, labelnew)
% =========================================================================
% Determine adaptively the best threshold to extract features from Y, 
% based on some sort of feature ranking score, using either:
% 1) the SVM, or alternatively the kNN algorithm (latter may be the better 
%    choice), 
% 2) one of the following performance cirteria:
%    a) the LOO crossvalidation performance on the training data, or 
%       alternatively ...
%    b) the criterion (Training Perf + Test Perf) / 2
% 
% Inputs:
% -------
% Y / Ynew              = [p * m] (preprocessed) training / test data
% labels / labelnew     = [p * 1] respective supervision info
% V                     = [m * 1] ranking list
% adap                  = algorithm params
%
% Outputs:
% --------
% adapthresh            = Optimal threshold of Y
% Crit                  = percentile
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 12/2013
function [adapthresh, Crit, mx, PredYnew, trainperf, perf] = nk_AdaptiveThreshold2(Y, labels, V, adap, Ynew, labelnew, sigma)

global FEATSEL MODEFL EVALFUNC SVM

V=abs(V);
nclass = size(V,2);
percvec = adap.perc_vec;
lpercvec = length(percvec);
perf = nan(lpercvec,nclass);
trainperf = nan(lpercvec,nclass);
sumfeat = zeros(lpercvec,nclass);
Crit = zeros(nclass,1);
if FEATSEL.showres
    figure(10);hold on
end
if size(V,2)>1, adapthresh=zeros(1,nclass); end
if size(labels,1) ~= size(Y,1), labels = labels'; end
colp = {'bx','rx','gx','mx','cx'};
colpt = {'b*','r*','g*','m*','c*'};
fprintf('\nAdaptively threshold feature volume:')

switch MODEFL
    case 'classification'
        func='max';
        if FEATSEL.BINMOD 
            if max(labels) == 2, labels(labels==2) = -1; end
            if max(labelnew) == 2, labelnew(labelnew==2) = -1; end
        end
    case 'regression'
        switch SVM.GridParam
            case {9,11,12,18}
                func='min';
            otherwise
                func='max';
        end
end

if isfield(adap,'SlackParam')
    SlackParam = num2str(adap.SlackParam,'%1.10f');
else
    SlackParam = '0';
end
if isfield(adap,'KernParam')
    KernParam = num2str(adap.KernParam,'%1.10f');
else
    KernParam = '0';
end

YnewPred = zeros(numel(labelnew), lpercvec, nclass);
YPred = zeros(numel(labels), lpercvec, nclass);

for curclass=1:nclass

    for k=lpercvec:-1:1

        % Apply current percentile threshold to volume
        if size(V,2) > 1, tV = V(:,curclass); else tV = V; end
        tV(isnan(tV)) = 0;
        thresh = percentile(tV,percvec(k));
        indthresh = tV>=thresh;
        tY = Y(:, indthresh);
        tYnew = Ynew(:, indthresh);
        sumfeat(k,curclass) = size(tY,2);
        %if k<lpercvec %%&& sumfeat(k+1,curclass)==sumfeat(k,curclass), 
        %    break, 
        %end
        
        % Train classifier with suprathreshold data
        switch adap.prog
            case 'LIBSVM'
                param.KernParam = KernParam;
                param.SlackParam = SlackParam;
            case 'LIBLIN'
                param.SlackParam = SlackParam;
            case 'IMRELF'
                param.distance = FEATSEL.imrelief.distance;
                param.sigma = sigma;
                param.V = tV(indthresh);
        end
        
        switch adap.method
            
            case 0 % (CV1-Training Perf. + CV1-Test Perf.) / 2 Criterion
                
                [YnewPred(:,k,curclass), YPred(:,k,curclass)] = eval_features(adap, tY, labels, tYnew, labelnew, param);
                trainval.acc = feval(EVALFUNC, labels, YPred(:,k,curclass));
                valt = feval(EVALFUNC, labelnew, YnewPred(:,k,curclass));
                fprintf('(CV1-Tr+CV1-Ts)/2: percentile = %1.3f%%, threshold=%1.3f, #features = %g, Train-Perf = %1.2f, CV1-Perf = %1.2f (%s)', ...
                    percvec(k), thresh, sumfeat(k,curclass), trainval.acc, valt(1), EVALFUNC);
                
            case 1 % LOO validation Criterion
                
                for q=1:lx
                    lx = size(tY,1); pred = zeros(lx,1);
                    ind = true(lx,1); ind(q) = false;
                    ltY = tY(ind,:); ltYnew = tY(q,:); lLabels = labels(ind);                    
                    pred(q) = eval_features(adap, ltY, lLabels, ltYnew, 0, param);
                end
                YPred(:,k,curclass) = pred;
                trainval.acc = feval(EVALFUNC, labels, pred); valt = trainval.acc;
                fprintf('LOO: percentile = %1.3f%%, threshold=%1.3f, #features = %g, LOO-Perf = %1.2f (%s)', ...
                    percvec(k), thresh, sumfeat(k,curclass), trainval.acc, EVALFUNC);
                
            case 2 % CV1-Test Perf. Criterion
                
                YnewPred(:,k,curclass) = eval_features(adap, tY, labels, tYnew, labelnew, param); 
                valt = feval(EVALFUNC, labelnew, YnewPred(:,k,curclass));
                if isnan(valt)
                    fprintf('\nNaN evaluation performance points to problems in CV structure setup.')
                    fprintf('\nProbably the CV1 test fold number is too high.')
                    error('Reset CV structure according to your group sizes.')
                end
                trainval.acc = valt(1);
                fprintf('CV1-Ts: percentile = %1.3f%%, threshold=%1.3f, #features = %g, CV1-Perf = %1.2f (%s)', ...
                    percvec(k), thresh, sumfeat(k,curclass), valt(1), EVALFUNC);
        end
        
        perf(k,curclass) = valt(1);
        trainperf(k,curclass) = trainval.acc;

    end
    
    if size(V,2) > 1,
        [mx1,in1] = feval(func,trainperf(:,curclass));
        [mx2,in2] = feval(func,perf(:,curclass));
        [mx,in] = feval(func,(trainperf(:,curclass)+perf(:,curclass))/2);
        adapthresh(curclass) = percentile(tV,percvec(in)); 
        Crit(curclass) = percvec(in);
        if FEATSEL.showres
            plot(sumfeat(in1,curclass), trainperf(in1,curclass),colp{curclass});
            plot(sumfeat(in2,curclass), perf(in2,curclass),colpt{curclass});
            ylim([20 100]);
            title(['CV1 performance of SVM problem #' num2str(curclass) ' depending on # features'])
            ylabel('performance')
            xlabel('# of features')
            drawnow
        end
    end
end

if size(V,2) == 1
    mperf = mean(perf,2);
    mtrainperf = mean(trainperf,2);
    [mx1,in1] = feval(func,trainperf);
    [mx2,in2] = feval(func,mperf);
    msumfeat = mean(sumfeat,2);
    [mx,in] = feval(func,(mperf+mtrainperf)/2);
    adapthresh = percentile(tV,percvec(in));
    Crit(1) = percvec(in);
    if FEATSEL.showres
        plot(msumfeat(in1), mtrainperf(in1),'kx');
        plot(msumfeat(in2), mperf(in2),'k*');
        ylim([20 100]);
        title('Mean CV1 performance of all SVM classifiers depending on # features')
        ylabel('performance')
        xlabel('# of features')
        drawnow
    end
end
if exist('YnewPred','var')
    PredYnew=zeros(size(YnewPred,1),nclass);
    for curclass=1:nclass
        PredYnew(:,curclass) = YnewPred(:,in,curclass);
    end
else
    PredYnew=[];
end
return

% _________________________________________________
function [Pnew, P] = eval_features(adap, Y, L, Ynew, Lnew, param)

switch adap.prog
                    
    case 'LIBSVM'

        fprintf('\nLIBSVM->')
        [~,model]   = nk_GetParam_LIBSVM(Y, L, param.SlackParam, param.KernParam, 1);
        Pnew        = nk_GetTestPerf_LIBSVM([],Ynew, Lnew, model);
        if nargout >1,
            P       = nk_GetTestPerf_LIBSVM([],Y,L, model);
        end

    case 'LIBLIN'

        fprintf('\nLIBLINEAR->')
        [~,model]   = nk_GetParam_LIBLIN( Y, L, param.SlackParam, [], 1);
        Pnew        = nk_GetTestPerf_LIBLIN([], Ynew, Lnew, model);
        if nargout >1,
            P       = nk_GetTestPerf_LIBLIN([], Y, L, model);
        end
        
    case 'kNNMEX'

        fprintf('\nkNN->')
        [~,model] = nk_GetParam_kNNMEX(tY, labels, [], [], 1);
        Pnew        = nk_GetTestPerf_kNNMEX([],Ynew, Lnew, model);
        if nargout >1 
            P       = nk_GetTestPerf_kNNMEX([],Y, L, model);
        end
        
    case 'IMRELF' 

        fprintf('\nIMRelief->')
        Pnew = IMRelief_Sigmoid_FastImple_Predict2(Y', L, Ynew', param.itV, param.distance, param.sigma);
        if nargout > 1
            P = IMRelief_Sigmoid_FastImple_Predict2(Y', L, Y', param.itV, param.distance, param.sigma);
        end
        P = P'; Pnew = Pnew';
end


return


