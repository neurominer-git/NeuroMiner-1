function [rs, ds] = nk_GetTestPerf_DECTRE(~, tXtest, ~, md, ~, ~)

[rs, ds] = md.predict(tXtest); 
ds = ds(:,2)-0.5;
% ind = tXtest(:,2)>md.binthr;
% rs = zeros(size(tXtest,1),1); 
% rs(ind)=1; rs(~ind)=-1; 
% ds = rs;

