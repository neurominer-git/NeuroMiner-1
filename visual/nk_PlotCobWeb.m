function [h, n_mat, ca, o] = nk_PlotCobWeb(mat, groupnames, ha)

cla(ha);
[M,N] = size(mat);
idx = eye(M,N);
n_mat =  (mat ./ sum(mat,2))*100;
mcl = n_mat(~idx);
Ax = cell(M,N); Ux = cell(M,N);
for i=1:M
    for j=1:N
        Ax{i,j} = sprintf('%s\n%s', groupnames{i},groupnames{j});
        Ux{i,j} = '';
    end
end
Ax = [Ax(~idx) Ux(~idx)];
% Plot cobweb graph
[h, ca, o] = spider([ones(numel(mcl),1)*100/N mcl ],'Misclassification web',[ zeros(numel(mcl),1) 50*ones(numel(mcl),1) ],Ax,[], ha);
o(1).LineStyle=':'; o(2).LineWidth=2.5;
aref = polyarea(o(1).XData,o(1).YData);
a2 = polyarea(o(2).XData,o(2).YData);
ca.Title.String=sprintf('[ RAR: %1.3f (%1.1f%% REF) ]',aref/a2, a2*100/aref);
ca.Title.Position = [0 1.26];
ca.Title.Visible='on';
