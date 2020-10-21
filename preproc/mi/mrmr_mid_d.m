function [fea] = mrmr_mid_d(d, f, K, bdisp)
% function [fea] = mrmr_mid_d(d, f, K)
% 
% MID scheme according to MRMR
%
% By Hanchuan Peng
% April 16, 2003
%
if ~exist('bdisp','var'), bdisp = 0; end;

nd = size(d,2);

t = zeros(nd,1);
for i=1:nd, 
   t(i) = mutualinfo(d(:,i), f);
end; 
%fprintf('calculate the marginal dmi costs %5.1fs.\n', cputime-t1);

[~, idxs] = sort(-t);
%fea_base = idxs(1:K);
fea = zeros(1, K);
fea(1) = idxs(1);

KMAX = min(1000,nd); %500

idxleft = idxs(2:KMAX);

k=1;
if bdisp==1,
fprintf('k=1 cost_time=(N/A) cur_fea=%d #left_cand=%d\n', ...
      fea(k), length(idxleft));
end;
curlastfea = 1;
for k=2:K,
   t1=cputime;
   ncand = length(idxleft);
   t_mi = zeros(ncand,1);
   c_mi = zeros(ncand,1);
   mi_array = zeros(ncand,curlastfea);
   for i=1:ncand,
      t_mi(i) = mutualinfo(d(:,idxleft(i)), f); 
      mi_array(idxleft(i),curlastfea) = getmultimi(d(:,fea(curlastfea)), d(:,idxleft(i)));
      c_mi(i) = mean(mi_array(idxleft(i), :)); 
   end;
    
   [~, fea(k)] = max(t_mi(1:ncand) - c_mi(1:ncand));

   tmpidx = fea(k); fea(k) = idxleft(tmpidx); idxleft(tmpidx) = [];

   if bdisp==1,
   fprintf('k=%d cost_time=%5.4f cur_fea=%d #left_cand=%d\n', ...
      k, cputime-t1, fea(k), length(idxleft));
   end;
   curlastfea = curlastfea + 1;
end;

return;

%===================================== 
function c = getmultimi(da, dt) 
for i=1:size(da,2), 
   c(i) = mutualinfo(da(:,i), dt);
end; 
    
