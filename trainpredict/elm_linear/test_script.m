function Params = test_script(X,Y, nHidden, C)
global VERBOSE
VERBOSE = false;
n = size(X,1); 
L = Y; L(L==2)=-1;
l=1;
for j=1:numel(C)
    for i=1:numel(nHidden)
        fprintf('\nnHidden = %g neurons; C = %1.3f: ', nHidden(i), C(j))
        kParams = zeros(10,1);
        for k = 1:10
            fprintf('.')
            trind = randperm(n,round(n/100*80));
            ind = false(n,1);
            ind(trind) = true;
            IN.ZeroOne = 2; IN.zerooutflag=1;
            tX1=X(ind,:); tX2=X(~ind,:);
            [tX1, IN] = nk_PerfScaleObj(tX1, IN);
            tX2 = nk_PerfScaleObj(tX2,IN); 
            TrL = Y(ind); TsL = L(~ind);
            tX1(isnan(tX1))=0;
            tX2(isnan(tX2))=0;
            %MATLAB SVM Implementation
            svmStruct = fitcsvm(tX1,TrL,'BoxConstraint',C(j),'KernelFunction','RBF','KernelScale','auto','OutlierFraction',0.05,'Solver','ISDA');
            T = predict(svmStruct,tX2);
            T(T==2)=-1;
            kParams(k) = BAC(TsL, T); 
            %LIBSVM
%             model = svmtrain312([],TrL,tX1,sprintf('-c %g -s 0 -t 0',C(j)));
%             T = svmpredict312(TsL,tX2,model);
%             T(T==2)=-1;
%             kParams(k) = BAC(TsL, T); 
            %LIBLIN
%             model = train_liblin(TrL,sparse(tX1),sprintf('-c %g -s 6',C(j)));
%             T = predict_liblin(TsL,sparse(tX2),model);
%             T(T==2)=-1;
%             kParams(k) = BAC(TsL, T); 
            %MEXELM
%             tX1=tX1'; tX2=tX2';
%             [inW, bias, outW] = mexElmTrain( tX1, TrL, nHidden(i), C(j));
%             elm_scores = mexElmPredict( inW, bias, outW, tX2 );
%             [~,T] = max(elm_scores); 
%             T(T==2)=-1;
%             kParams(k) = BAC(TsL, T'); 
            clear IN
        end
        Params(l) = nm_nanmean(kParams);
        fprintf('\nBAC = %g',Params(l))
        l=l+1;
    end
end
figure;plot(Params)