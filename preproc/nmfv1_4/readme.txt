+------------------------------------------------+
| Written by: Vijay Kumar                        |
|             Queen Mary, University of London   |
|                                                |
+------------------------------------------------+
This code implements the paper :

Vijay. K, I. Kotsia and I. Patras, Max-margin Semi-NMF, in British Machine Vision Conference, 2011.

Input parameters:
1) X (m x n): training matrix where m is the dimension of the features and n is the number of training features.
2) r : reduced dimension
3) labl : label for the training set
4) C : SVM parameter [0.1 to 100]
5) K : scaling parameter [10 to 10000]

Output parameters:
1) G (m x r): learned bases
2) H (r x n): learned coefficients 
3) w (r x 1): svm hyperplane parameter

Requirements : 1) Requires libsvm package installed
                  http://www.csie.ntu.edu.tw/~cjlin/libsvm/
               2) Also requires package 'Quadratic programming in C' 
                  http://sigpromu.org/quadprog/

If you find this code useful, please cite it as
Vijay. K, I. Kotsia and I. Patras, Max-margin Semi-NMF, in British Machine Vision Conference, 2011.
