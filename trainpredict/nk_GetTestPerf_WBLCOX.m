function [rs, ds, md] = nk_GetTestPerf_WBLCOX(tXtest, md)
global SVM
lx = size(tXtest,1); 
if isfield(md,'interval')
    interval = md.interval;
else
    interval = mean(md.t0); % median time
end

ds = 1 - exp( - ((interval/md.rho) .^ md.nu) * exp(tXtest*md.beta) );

if isfield(md,'cutoff')
    md.cutoff.prob = percentile(ds,md.cutoff.val);
else
    md.cutoff.prob = percentile(ds,50); % median time
end

% Compute predicted outcomes at given outcomes but do not change
% probabilities
dds = nk_CalibrateProbabilities(ds, md.cutoff.prob ); rs = sign(dds);

% Compute predicted time-2-event
if nargout>2 
    if SVM.WBLCOX.predict_time
        md.predicted_time = arrayfun( @(i) wc_predict(tXtest(i,:), md), 1:lx )';
    else
        md.predicted_time = nan(size(tXtest,1),1);
    end
end