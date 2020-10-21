load faces.mat;
X = X / 255;

disp 'Press any key to show example faces in the input dataset.';
pause;

ShowNMFBasis(X(1:20:2000,:)', 32, 32, 2, true);

disp 'Press any key to run the PNMF algorithm based on the Euclidean distance.';
disp '   W = pnmfeu(X'', 16)'
pause;

W = pnmfeu(X', 16);

disp 'Learning is finished.'
disp 'The resulting matrix W is highly sparse'.
disp 'Press any key to show the resulting columns of W as basis images.';
pause;

ShowNMFBasis(W, 32, 32, 2, true);

disp ' '
disp 'Next we demonstrate application of PNMF on sparse data. Press any key to continue.'
pause;

clear;
load cisi.mat;

disp 'cisi data loaded.'
disp 'Press any key to spy the input cisi data matrix.'
disp 'Blue indicates non-zeros.'
pause;

figure; spy(X);

disp 'Press any key to run PNMF algorithm based on generalized KL-divergence (I-divergence).'
disp '    W = pnmfkl_sparse(X, 10)'
pause;

W = pnmfkl_sparse(X, 10);

disp 'Learning is finished.'
disp 'The resulting matrix W is highly sparse, which can be seen by zeroing out the very small entries.'
disp 'W(W<1e-20) = 0'
disp 'Press any key to display the matrix W.'
pause

W



