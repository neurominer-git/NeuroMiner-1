function [ts, ds, rs, hrs, hds, hrx, hdx, hts_rs, hts_ds] = nk_GetTestParam(Y, nclass, model, features, weights)
%
% function [ts, ds, rs, hrs, hds, hrx, hdx, hts_rs, hts_ds] =
% nk_GetTestParam(Xtest, Ytest, model, features)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% INPUTS
% ------
% Y:        Input feature structure
% Ytest:    test data values to predict
% model:    learned prediction rule
% features: subspaces to be used from Xtest

% OUTPUTS
% -------
% ts:       Test data performance (see SVM.GridParam)
% ds:       Decision values (only in SVM)
% rs:       Predicted target values
% hrs:      Aggregated predicted target values (ensemble learning)
% hds:      Aggregated decision values (ensemble learning, SVM only)
% hrx:      Ensemble decision using majority voting: 
%           sign(sum(H)) (for brinary classification)
% hdx:      Ensemble decision using decision value averaging (SVM only)
% hts_rs:   Ensemble performance using majority voting
% hrs_ds:   Ensemble performance using decision value averaging
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 12/2009

global MODEFL SVM RFE MULTI

[ic, jc] = size(model); [np, nf, nv] = size(Y.Ts);

if jc>1 && ic>1
    
    rs = cell(ic,jc);
    ds = rs;
    ts = rs;
 
    for i=1:ic
        
        for j=1:jc
        
            tX = cell(1,nv);
            tXtest = tX;
            
            % 1) Compute binary classifier / predictor peformance on CV2 test data
            for curclass=1:nclass
                for v=1:nv %loop through variates                        
                    if ~iscell(Y.Tr{i,j,v})
                        XTr = Y.Tr{i,j,v}(Y.TrInd{i,j}{curclass},:);
                        XCV = Y.CV{i,j,v}(Y.CVInd{i,j}{curclass},:);
                        XTest = Y.Ts{i,j,v}(Y.TsInd{curclass},:);
                    else
                        XTr = Y.Tr{i,j,v}{curclass}(Y.TrInd{i,j}{curclass},:);
                        XCV = Y.CV{i,j,v}{curclass}(Y.CVInd{i,j}{curclass},:);
                        XTest = Y.Ts{i,j,v}{curclass}(Y.TsInd{curclass},:);
                    end
                    if RFE.ClassRetrain, 
                        tX{v} = [XTr; XCV]; 
                    else
                        tX{v} = XTr; 
                    end
                    tXtest{v} = XTest;
                end
                Ytest = Y.TsL{curclass};
                
                % This for binary optimized features subspaces
                if ~MULTI.flag || (MULTI.flag && ~MULTI.train)
                    F = squeeze(features(i,j,curclass))';
                else
                    F = squeeze(features(i,j))';
                end
                
                [ts{i,j,curclass}, rs{i,j,curclass}, ds{i,j,curclass}] = nk_GetTestPerf(tXtest, Ytest, F, model{i,j}, tX);
            end
        end
    end
else
    
    if ~exist('features','var'), F = []; else F = features; end
    if ~isempty(XCV) && RFE.ClassRetrain, X = [XTr; XCV]; else X = XTr; end
    [m,n] = size(Xtest);
    [ts, rs, ds] = get_test_perf(Xtest, Ytest, F, model, X);

end

%% Compute classifier / predictor performance on (hold-out) data
if iscell(ts)
    
    ENSMETHOD = 1;
    [hts, hrs, hds] = ensemble_classout(ts,rs,ds, weights);
    
    switch MODEFL
        
        case 'classification'
            
            switch ENSMETHOD
                
                case 1 % This is the simple majority vote
                    
                    % Hard label decision:
                    hrx = sign(sum(hrs,2));
                    
                    % Soft decision value / probability decision:
                    if SVM.RVMflag || strcmp(SVM.Optimization.b,' -b 1')
                        mhds = mean(hds,2); 
                        hdx = zeros(size(mhds));
                        hdx(mhds > 0.5) = 1;
                        hdx(mhds < 0.5) = -1;
                    else
                        hdx = mean(hds,2);
                    end
                    
                case 2 % Product majority vote
                    hrx = sign(prod(hrs,2));
                    hdx = prod(hds,2);
                
                case 3 % Error Correcting Output Codes
                    coding=1; decoding=1;
                    classes = ones(1,size(hrs,2));
                    %hrs(hrs==-1) = ;
                    hrx = nk_ErrorCorrOutCodes(hrs, classes, coding, decoding);
                    hrx(hrx==2) = -1; hdx = hrx;
            end
               
            hts_rs  = get_test_param(SVM.GridParam, Ytest, hrx);
            hts_ds  = get_test_param(SVM.GridParam, Ytest, sign(hdx));
            
        case 'regression'
            hrx     = mean(hrs,2);
            hdx     = hrx;
            hts_rs  = get_test_param(SVM.GridParam, Ytest, hrx);
            hts_ds  = hts_rs;
    end
    
end

return

% _______________________________________
function [hts, hrs, hds] = ensemble_classout(ts, rs, ds, w)
    
[ix,jx] = size(rs);cnt=1; hts = []; hrs = []; hds = []; hdw = [];
for i=1:ix
    for j=1:jx  
        hts = [hts; ts{i,j}];
        hrs = [hrs rs{i,j}];
        hds = [hds ds{i,j}];
        hdw = [hdw w{i,j}'];
        cnt=cnt+1;
    end
end
m = size(hrs,1);
wx = repmat(hdw,m,1);
hrs = hrs.*wx;
hds = hds.*wx;
return


