function Xts = nk_PLSpredict(brainlv_tr,Ytr,Xtr,Yts)

mx = max(Xtr);
Xd = single(zeros(size(Xtr,1),mx));
dx = Xtr(1); cnt=1;Xd(1)=dx;
for i=2:length(Xtr)
    if Xtr(i) > dx,
        dx=Xtr(i);
        cnt=cnt+1;
    end
    Xd(i,cnt)=1;
end
%[mtr,ntr]=size(Ytr);
mts=size(Yts,1);
% Mean centering
mYtr = mean(Ytr);
mXtr = mean(Xd);

%Ytr_mce = Ytr-(ones(mtr,1)*mYtr);
Yts_mce = Yts-(ones(mts,1)*mYtr);

Xts=Yts_mce*brainlv_tr;				              
Xts=Xts+mXtr;

return