function [ R, P ] = nk_SimNfeatNcaseNmarkers(Nfeats, Nmarkers, Ncases, eventprob, noiseiter)

P = allcomb2(Nfeats, Nmarkers, Ncases, eventprob, noiseiter);
nP = size(P,1);
R = zeros(nP,1);
return
figure;
poolobj=parpool(5);
for i=1:nP
    fprintf('\nWorking on parameter combination %g/%g: N_feats = %g, N_markers = %g, N_cases = %g, eventrate = %g, noise iterations = %g', i, nP, P(i,:));
    R(i) = compute_perf(P(i,1), P(i,2), P(i,3), P(i,4), P(i,5));
    plot(1:i,R(1:i),'b-'); drawnow
end
delete(poolobj);

function R = compute_perf(nf, mr, nc, er, nnois)
global SVM
nm = mr;
nc1 = ceil(nc*er);
nc2 = nc-nc1;
% Create label
L = [ones(nc1,1); -1*ones(nc2,1)];
% Create marker matrix
Mx = zeros(nc,nm);
for i=1:nm
    Mx(:,i) = [rand(nc1,1); -1*rand(nc2,1)];
end
M = [Mx rand(nc,nf-nm)];

RAND.OuterPerm = 1;
RAND.InnerPerm = 1;
RAND.OuterFold = 5;
RAND.InnerFold = 5;
RAND.Decompose = 1;
SVM.prog = 'LIBSVM';
SVM.(SVM.prog).Weighting = 1;
SVM.(SVM.prog).Optimization.b = 0;
MODEFL = 'classification';
LIBSVMTRAIN = 'svmtrain312';
LIBSVMPREDICT = 'svmpredict312';
CMDSTR.WeightFact = 1;
Lx = L; Lx(L==-1)=2;
cv = nk_MakeCrossFolds(Lx, RAND, 'classification', [], {'A','B'} );
[pCV,nCV] = size(cv.TrainInd);
reps = 10;
R = zeros(reps,pCV*nCV);
for k=1:reps
    Ystar = M;
    for i=1:nnois
        % Create random matrix
        Y = rand(nc, nf);
        % Marker matrix with random data
        Ystar = Ystar.*Y;
    end
    Ystar = nk_PerfScaleObj(Ystar);
    Ri = zeros(pCV,nCV);

    for j=1:pCV
        parfor i=1:nCV
            Tr = cv.TrainInd{j,i};
            Ts = cv.TestInd{j,i};
            cmd = '-s 0 -t 0 -c 1';
            cmd = nk_SetWeightStr(SVM, MODEFL, CMDSTR, L(Tr), cmd);
            model = feval( LIBSVMTRAIN, [], L(Tr), Ystar(Tr,:), cmd);
            [~, ~, ds] = feval( LIBSVMPREDICT, L(Ts), Ystar(Ts,:), model, sprintf(' -b %g',SVM.LIBSVM.Optimization.b)); 
            Ri(j,i) = BAC(L(Ts),ds);
        end
    end
    R(k,:) = Ri(:);
end
R = nm_nanmean(R(:));