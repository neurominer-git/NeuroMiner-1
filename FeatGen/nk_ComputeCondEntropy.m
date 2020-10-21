function vargout = nk_ComputeCondEntropy(dat, ind, kstr, gstr)

global SVM

[mapY,mapYnew] = nk_Preprocess(dat, 50, ind);
[nperms, fold] = size(dat.cv.TrainInd);
if isstruct(mapY)
    dims = size(mapY.Y{1,1},2);
else
    dims = size(mapY,2);
end
perfTrain = zeros(dims,nperms*fold);
perfCV = zeros(dims,nperms*fold);
perfTest = zeros(dims,nperms*fold);

for h=1:dims
    l=1;
    fprintf('\nDim: %g',h)
    for i=1:nperms
        
        for j=1:fold
            
            if isstruct(mapY)
                tY = mapY.Y{i,j}(dat.cv.TrainInd{i,j},1:h);
                tYcv = mapY.Y{i,j}(dat.cv.TestInd{i,j},:);
                tYnew = mapYnew{i,j};
            else
%                 fprintf('\nInd-Max: %g',max(dat.cv.TrainInd{i,j}));
%                 fprintf('\nInd-Min: %g',min(dat.cv.TrainInd{i,j}));
%                 fprintf('\nDim: %g', h)
%                 fprintf('\nY-x: %g, Y-y: %g)',size(mapY,1),size(mapY,2))
                tY = mapY(dat.cv.TrainInd{i,j},1:h);
                tYcv = mapY(dat.cv.TestInd{i,j},:);
                tYnew = mapYnew;
            end
            
            lY = dat.label(dat.cv.TrainInd{i,j});
            lYcv = dat.label(dat.cv.TestInd{i,j});
            lYnew = dat.labelnew;
            
            %[RF, model] = nk_RFEfold(tY, lY, tYcv, lYcv, SVM, kstr, gstr, 0, ones(1:h,1));
                        
            % Build SVM using symbolized features
            [Train, model] = get_param(tY, lY, SVM, kstr, gstr,1);
            perfTrain(h,l)=Train.val;%RF.AchievedTrainParameter;

            % Predict class labels of CV & test data
            perfCV(h,l) = gettestvals(SVM, tYcv, lYcv, model, 1:h);

            
            % Predict class labels of CV & test data
            perfTest(h,l) = gettestvals(SVM, tYnew, lYnew, model, 1:h);
            l=l+1;
        end
        
    end
end

vargout.Perf_Train = mean(perfTrain,2);
vargout.Perf_Test = mean(perfTest,2);
vargout.Perf_CV = mean(perfCV,2);
figure(3)
%plot(2,1,2);
plot(vargout.Perf_Train,'b-'), hold on
plot(vargout.Perf_Test,'b:') 
plot(vargout.Perf_CV,'b--')
return