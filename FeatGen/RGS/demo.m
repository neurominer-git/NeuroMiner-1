function demo(func, k);

% demo;
% demo(func);
% demo(func, k);
%
% Demo for using the algorithms in this directory, developed or used in
%
% A. Navot, L. Shpigelman, N. Tishby, E. Vaadia. Nearest Neighbor Based Feature Selection for Regression and its
% Application to Neural Activity. Submitted to NIPS 2005.
%
% the demo demonstrats the effect of feature selection and shows the
% results of the various selection algorithms on the given 2-dimesional
% function
%
% input: func: the target function to be estimated. valid values are
%              'sin_x1' (default), 'sin_x1_2', 'square_x1'
%        k: number of neighbors. (default is 3)
%
% 
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Written by Amir Navot & Lavi Shpigelman               
%% Date: June 3, 2005                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <2,
    k = 4;
end

if nargin <1,
    func = 'sin_x1';
end

X_train = rand(100,2);

Ytr = feval(func, X_train(:,1), X_train(:,2));
Y_train = Ytr(:) + rand(100,1)/5;

fh = figure; 
ss = get(0,'ScreenSize');
set(gcf, 'Position', ss/1.5 + [0 ss(4)/3 0 0] )

x1 = 0:0.01:1;
x2 = x1;
[X1, X2] = meshgrid(x1, x2);

X_test = [X1(:), X2(:)];
Y = feval(func, X1, X2);
Y_test = Y(:);

%figure; 
subplot(221)
myplot(X1, X2, Y);

title('The Function We Try to Estimate');

subplot(222)
plot3(X_train(:,1), X_train(:,2), Y_train, 'b.');

title('The Training Set');
grid on


extra_param.k = k;
estY = kNNreg(X_test, X_train, Y_train, [1 1], extra_param);
estYr = reshape(estY, length(x1), length(x1));


subplot(223);
myplot(X1,X2,estYr);

title(['Estimated Function, using both features'])
xlabel(['mean square error = ' num2str(mean((estY-Y_test).^2))]);


% feature selection:
%RGS
extra_param.num_starts = 1;
extra_param.epochs = 1;
weights = RGS(X_train, Y_train, extra_param);
[w ord] = sort(-weights);
disp(['The weights that RGS assigned are: (' num2str(weights(1)) ', ' num2str(weights(2)) ')']);

% SKS
scores = SKS(X_train, Y_train, extra_param);
[w ord] = sort(-scores);
disp(['The scores that SKS assigned are: (' num2str(scores(1)) ', ' num2str(scores(2)) ')']);

% correletion coeffitient
scores = corrcoef_select(X_train, Y_train);
[w ord] = sort(-scores);
disp(['The scores that Corr_Coef assigned are: (' num2str(scores(1)) ', ' num2str(scores(2)) ')']);


% correletion coeffitient
scores = infoGain(X_train, Y_train);
[w ord] = sort(-scores);
disp(['The scores that infoGain assigned are: (' num2str(scores(1)) ', ' num2str(scores(2)) ')']);

% forword selection
selected_features = fwd_sel(X_train, Y_train, extra_param);
disp(['The features selected by  fwd_sel are: ' num2str(selected_features)]);

%first feature only
w = [1 0];
estY = kNNreg(X_test, X_train, Y_train, w , extra_param);
estYr = reshape(estY, length(x1), length(x1));

subplot(224);
myplot(X1,X2,estYr);
title(['Estimated Function, using first feature only'])
xlabel(['mean square error = ' num2str(mean((estY-Y_test).^2))]);
%disp(['kNNreg using the best feature only. MSE = ' num2str(mean((estY-Y_test).^2)) ' (hit Enter)']);


function Y = square_x1(X1, X2)

Y = X1.^2;


function Y = sin_x1(X1, X2)

f = 4*pi;
Y = sin(f*X1+pi/2);

function Y = sin_x1_2(X1, X2)

f = 2*pi;
Y = sin(f*X1+pi/2);

function myplot(X1,X2,Y)
surf(X1,X2,Y); shading flat