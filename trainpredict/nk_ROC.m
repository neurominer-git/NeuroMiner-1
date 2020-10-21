function r = nk_ROC(labels,scores,posclass)
[r.xroc, r.yroc, r.t, r.auc, r.optpt] = perfcurve(labels,scores,posclass);
end