This is the Sparse Representation Toolbox in MATLAB version 1.9.

Sparse Reprsentation:
Kernel sparse coding: KSRSC
Kernel dictionary learning: KSRDL
Versatile sparse matrix factorization (VSMF): vsmf 

Classifiers:
The classification methods include
NNLS: nnlsClassifier, bootstrapnnlsClassifier
SRC1: src
SRC2: SRC2
MSRC: computeMetaSample and SRC2
MNNLS(Meta-sample based NNLS): computeMetaSample and nnlsClassifier
Least regression classifier (LRC): lrc
Nearest centroid: nearestCentroidTrain, nearestCentroidPredict
C-SVM and nu-SVM: svmTrain, softSVMTrain2, softSVMPredict2
High dimensional linear machine: khdlmTrain, khdlmPredict
Unified interface: classificationTrain, classificationPredict

Optimization:
NNQP using active-set algorithm: NNQPActiveSet (l1NNQPActiveSet)
l1QP using active-set algorithm: l1QPActiveSet
NNQP using MSO algorithm: NNQPMSO, NNQPSMOMulti
l1QP using MSO algorithm: l1QPMSO, l1QPSMOMulti
NNQP using interior-point algorithm: NNQPIP, NNQPIPMulti
l1QP using interior-point algorithm: l1QPIP, l1QPIPMulti
l1QP using proximal algorithm: l1QPProximal, l1QPProximalMulti

Statistical Test:
Find significant accuracy/error rate: significantAcc
Fitting learning curves: learnCurve
Friedman test coupled with Nemenyi test: FriedmanTest
Plot the result of Nemenyi test: plotNemenyiTest
Leave-M-Out: leaveMOut

Others:
Plot the accuracy / running time of multiple classifiers on multiple datasets: plotBarError.m

Requirement:
1. The codes are written in MATLAB 7.11(2010b). It may also work in earlier version, for example MATLAB 7.9(2009a).
2. Some functionalities need my NMF toolbox: https://sites.google.com/site/nmftool/.
3. For using SRC1: you have to install L1_magic which is available at http://users.ece.gatech.edu/~justin/l1magic.
4. For using SRC2: you have to install l1_ls which is available at http://www.stanford.edu/~body/l1_ls, if you choose the 'interiorPoint' method.

Usage:
Please see the files example*.m for how to use them.
For more parameter seting, please use command HELP in the command line, for example help nnlsClassifier. 

Documentation:
There is no complete documentation by now, I am still working on this. The following references should be very useful for understanding this toolbox.

References:
[1] Li Y, and Ngom A. Versatile sparse matrix factorization: theory and applications Neurocomputing 145:23-29. (2014)
[2] Li Y, Caron R.J., Ngom A. A decomposition method for large-scale sparse coding in representation learning. IJCNN, 2014, pp. 3732-2738
[3] Li Y, Ngom A. Sparse representation approaches for the classification of high-dimensional biological data. BMC Syst Biol 7 Suppl 4:S6. (2013) PMID 24565287
[4] Li Y, Ngom A. The non-negative matrix factorization toolbox for biological data mining. Source Code Biol Med 8(1):10. (2013) PMID 23591137
[5] Li Y, Ngom A. Non-Negative Least Squares Methods for the Classification of High Dimensional Biological Data. IEEE/ACM Trans Comput Biol Bioinform (2013) PMID 23568964
[6] Li Y, Ngom A. Classification approach based on non-negative least squares Neurocomputing 118:41-57. (2013)
[7] Li Y, Ngom A. Supervised Dictionary Learning via Non-negative Matrix Factorization for Classification," ICMLA 2012, pp. 439-443.

Contact:
Should you have any question and problem about this toolbox, please contact 
Yifeng Li
CMMT, UBC, Vancouver, Canada
yifeng@cmmt.ubc.ca
li11112c@uwindsor.ca
yifeng.li.cn@gmail.com

March 02, 2015
