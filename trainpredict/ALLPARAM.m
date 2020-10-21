% =========================================================================
% FORMAT param = ALLPARAM(expected, predicted, invflag, tiedfun)
% =========================================================================
% Compute all classification parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 10/2017

function param = ALLPARAM(expected, predicted, invflag, tiedfun )
if isempty(expected), param = []; return; end
if ~exist('invflag','var') || isempty(invflag), invflag = false; end
ind0 = expected ~=0 & ~isnan(predicted);
expected = expected(ind0); 
predicted = predicted(ind0);

if any(~predicted)
    mess = sprintf('%g tied predictions dectected.',sum(~predicted));
    if exist('tiedfun','var')
        switch tiedfun
            case 'neg'
                predicted(~predicted) = -1*realmin;
                mess = sprintf('%s Correcting toward -1*realmin',mess); 
            case 'pos'
                predicted(~predicted) = realmin; 
                mess = sprintf('%s Correcting toward +1*realmin',mess); 
        end
    else
        mess = sprintf('%s Skipping subjects with tied predictions in the computation.',mess); 
    end
    warning('%s',mess);
end

if invflag, predicted = -1*predicted; expected = -1*expected; end

TP = sum( predicted > 0 & expected > 0 );
FP = sum( predicted > 0 & expected < 0 );
TN = sum( predicted < 0 & expected < 0 );
FN = sum( predicted < 0 & expected > 0 );

param.TP = TP;
param.TN = TN;
param.FP = FP;
param.FN = FN;
param.P = TP + FP;
param.N = TN + FN;
param.acc   = ((TP + TN) / (TP+FP+TN+FN)) * 100 ;
param.sens  = (TP / (TP + FN)) * 100;
param.spec  = (TN / (TN + FP)) * 100;
param.FPR   = (FP / (TN + FP)) * 100;
param.PPV   = (TP / (TP + FP)) * 100;
param.NPV   = (TN / (TN + FN)) * 100;

sens = param.sens/100; spec = param.spec/100;

TPrate = TP / ( TP + FN); 
TNrate = TN / ( TN + FP); 
param.GMean = sqrt(TPrate + TNrate);

% Calculate further measures
if numel(unique(expected))>1 , param.AUC = fastAUC(expected, predicted, 1); end
param.BAC = (param.sens +  param.spec)/ 2;
param.Fscore = (2 * TP) / (2 * TP + FP + FN) ;
param.MCC = (TP * TN - FP * FN) / sqrt( (TP+FP) * (TP+FN) * (TN+FP) * (TN+FN) );
if param.spec == 100,
    % Avoid INF
    param.pLR = param.sens / (100 - 99.999999);
else
    param.pLR = param.sens / (100 - param.spec);
end
param.nLR = (100 - param.sens) / param.spec;
param.PSI = param.PPV + param.NPV - 100;
param.NNP = 1 / (param.PSI / 100);
param.NND = 1 / (sens - (1 - spec ));
param.Youden = ( sens + spec ) - 1;
param.DOR = sens/(1-spec)/((1-spec)/sens);
param = KAPPA(TP,TN,FP,FN,param);
end