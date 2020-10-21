function [fig, targets_means, predictions_means] = plotreliability(targets, predictions, show_fig)
%PLOTRELIABILITY Reliability diagram for calibration of two-class predictor.
%  PLOTRELIABILITY(TARGETS, PREDICTIONS) uses TARGETS and PREDICTIONS
%  vectors to compute TARGETS_MEANS and PREDICTIONS_MEANS. FIG is a
%  rudimentary reliability diagram plotted using these means.
%
%  Values of TARGETS are either 0 or 1, representing the two classes. 
%  Values of PREDICTIONS are probabilities of the positive class (1), 
%  previously obtained using a probabilistic predictive model. TARGETS and
%  PREDICTIONS are row or column vectors of the same size.
%
%  The predictions are discretized into ten bins. Cases with predicted
%  value between 0 and 0.1 fall in the first bin, between 0.1 and 0.2 in
%  the second bin, etc. For each bin, the mean predicted value
%  (PREDICTIONS_MEANS) and the true fraction of positive cases
%  (TARGETS_MEANS) are computed. These points are then plotted on the X and
%  Y axes respectively. If the model is well calibrated, the points will
%  fall near the diagonal line.
%
%  SHOW_FIG is an optional logical input argument. Its default value is
%  true. Setting it to false prevents the automatic display of the
%  reliability diagram.
%
%  Syntax:
%    plotreliability(targets, predictions);
%    [fig, targets_means, predictions_means] = plotreliability(targets, predictions, false);
%
%  Remarks:
%    This function does not allow the number of bins or the width of each 
%    bin to be changed.
%
%  References:
%    [1] Alexandru Niculescu-Mizil and Rich Caruana (2005) Predicting Good 
%        Probabilities With Supervised Learning, in Proceedings of the 22nd
%        International Conference on Machine Learning.
%        See section 4 (Qualitative Analysis of Predictions).
%
% Version: 20091103
% MATLAB version: 7.9.0.529 (R2009b)
%{
Keywords:
  reliability, reliability diagram, reliability plot, calibration, 
  machine learning
%}
%% Parse input
if nargin < 3, show_fig = true; end
%% Compute means
bins = max(ceil(predictions * 10), 1);
[predictions_means, targets_means] = deal(zeros([10, 1]));
for bin = 1:10
    predictions_in_bin = predictions(bins == bin);
    predictions_means(bin) = mean(predictions_in_bin);
    
    targets_in_bin = targets(bins == bin);
    targets_means(bin) = sum(targets_in_bin)/numel(targets);
end
%% Plot figure
fig = figure('Visible', 'off');
axes1 = axes('Parent', fig);
box(axes1, 'on')
hold(axes1, 'all')
plot(axes1, [0, 1], [0, 1], 'Color', 'black')
plot(axes1, predictions_means, targets_means, 'Color', 'blue', ...
     'Marker', 'o')
xlabel(axes1, 'Mean Predicted Value')
ylabel(axes1, 'Fraction of Positives')
title(axes1, 'Reliability Diagram')
hold(axes1, 'off')
if show_fig, set(fig, 'Visible', 'on'), end