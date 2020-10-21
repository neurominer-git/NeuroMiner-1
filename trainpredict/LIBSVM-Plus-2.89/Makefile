# This Makefile is used under Linux

#MATLABDIR = /home/nkoutsou/matlab2009B/
MATLABDIR = /usr/local/MATLAB/R2012b/
CXX = g++
#CXX = g++-4.1
CFLAGS = -Wall -O3 -fPIC -fopenmp -I $(MATLABDIR)/extern/include

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC=$(CXX) CXX=$(CXX) CFLAGS="$(CFLAGS)" CXXFLAGS="$(CFLAGS)"
# comment the following line if you use MATLAB 7.3 (or above) on 32-bit computer
MEX_OPTION += -largeArrayDims -lgomp
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)

all:	matlab

matlab:	binary

binary: svmpredict289PLUS.$(MEX_EXT) svmtrain289PLUS.$(MEX_EXT) read_sparse289PLUS.$(MEX_EXT) svmscale289PLUS.$(MEX_EXT)

svmpredict289PLUS.$(MEX_EXT):     svmpredict289PLUS.c svm.h svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmpredict289PLUS.c svm.o svm_model_matlab.o

svmtrain289PLUS.$(MEX_EXT):       svmtrain289PLUS.c svm.h svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmtrain289PLUS.c svm.o svm_model_matlab.o

read_sparse289PLUS.$(MEX_EXT):	read_sparse289PLUS.c
	$(MEX) $(MEX_OPTION) read_sparse289PLUS.c

svmscale289PLUS.$(MEX_EXT): svmscale289PLUS.c
	$(MEX) $(MEX_OPTION) svmscale289PLUS.c -o svmscale289PLUS

svm_model_matlab.o:     svm_model_matlab.c svm.h
	$(CXX) $(CFLAGS) -c svm_model_matlab.c

svm.o:  svm.cpp svm.h
	$(CXX) $(CFLAGS) -c svm.cpp

#test.$(MEX_EXT):	testgradient.c svm.h svm.o svm_model_matlab.o
#	$(MEX) $(MEX_OPTION) testgradient.c svm.o svm_model_matlab.o

clean:
	rm -f *~ *.o *.mex* *.obj
