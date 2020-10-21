function [condEntY, symY, binY, condEntYnew, symYnew] = nk_FilterCondEntropy(Y, Ynew, labels, labelsnew, SymFlag)

ix = size(Y,2);
condEntY = zeros(ix,1);

if nargin < 5
	SymFlag = 1;
end
	
if SymFlag
	% first symbolize features of training data
	[symY, binY, symYnew] = nk_Symbolize(Y, Ynew, 5, 5, 15);
else
	symY = Y;
	symYnew = Ynew;
end

if ~isempty(symYnew)
    condEntYnew = zeros(ix,1);
end

% Compute conditional entropy of each symbolized feature and the class
% labels for the training and test data

for i=1:ix

    fprintf('\nComputing conditional entropy of feature %g',i)
    if any(symY(:,i))
        condEntY(i) = condentropy(symY(:,i),labels);
        if ~isempty(symYnew)
            condEntYnew(i) = condentropy(symYnew(:,i),labelsnew);
        end
    else
        condEntY(i) = NaN;
    end
end


%mRMRTrain = mrmr_mid_d(symY,labels,ix);
[scondEntY, sI] = sort(condEntY,'ascend');
figure(2)
subplot(2,1,1); plot(condEntY,'b.') 
subplot(2,1,2); plot(scondEntY,'b.')

if strcmp(spm('ver'),'SPM5'), 
    spm5flag = 1; 
else
    spm5flag = 0;
end

if spm5flag
    Vm = spm_vol(spm_select(1,'image','select brainmask'));
else
    Vm = spm_vol(spm_get(1,'IMAGE','select brainmask'));
end

indvol=[];
% Read brainmask into vector format
for sl=1:Vm.dim(3)
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    %M1  = Vm.mat\V.mat\M;
    mask_slice = spm_slice_vol(Vm,M,Vm.dim(1:2),1);
    ind0 = find(mask_slice > 0.5);
    ind = ind0 + (sl - 1)*prod(Vm.dim(1:2));
    indvol = [indvol; ind];
    clear mask_slice
end

% Copy volume structure to SVVol
Condvol       = Vm;
Condvol.dim   = [Vm.dim(1:3) 4];
Condimg        = zeros(Vm.dim(1:3));

Condvol.fname = [pwd filesep 'CondVol.img'];
Condimg(indvol) = condEntY;
spm_write_vol(Condvol,Condimg);
return