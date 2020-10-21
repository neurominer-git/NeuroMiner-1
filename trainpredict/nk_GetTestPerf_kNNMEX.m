function [rs, ds] = nk_GetTestPerf_kNNMEX(~, tXtest, ~, md, ~, ~)

[target_test, predict_test] = fknn(md.X, md.Y, tXtest, [], md.kNN, 0);
if md.sfl, target_test(target_test==2) = -1; end
rs = target_test(:,1); ds = predict_test(:,1)-0.5;
            
end