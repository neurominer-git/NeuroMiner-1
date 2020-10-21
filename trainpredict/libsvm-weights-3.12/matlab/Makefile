# This Makefile is used under Linux

MATLABDIR = /usr/local/MATLAB/R2012b
# for Mac
# MATLABDIR ?= /opt/local/matlab

CXX ?= g++
#CXX = g++-4.1
CFLAGS = -Wall -Wconversion -O3 -fPIC -fopenmp -I$(MATLABDIR)/extern/include -I..

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)"
# comment the following line if you use MATLAB on 32-bit computer
MEX_OPTION += -largeArrayDims -lgomp
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)

OCTAVEDIR ?= /usr/include/octave
OCTAVE_MEX = env CC=$(CXX) mkoctfile
OCTAVE_MEX_OPTION = --mex
OCTAVE_MEX_EXT = mex
OCTAVE_CFLAGS = -Wall -O3 -fPIC -I$(OCTAVEDIR) -I..

all:	matlab

matlab:	binary

octave:
	@make MEX="$(OCTAVE_MEX)" MEX_OPTION="$(OCTAVE_MEX_OPTION)" \
	MEX_EXT="$(OCTAVE_MEX_EXT)" CFLAGS="$(OCTAVE_CFLAGS)" \
	binary

binary: svmpredict312.$(MEX_EXT) svmtrain312.$(MEX_EXT) libsvmread312.$(MEX_EXT) libsvmwrite312.$(MEX_EXT)

svmpredict312.$(MEX_EXT):     svmpredict312.c ../svm.h ../svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmpredict312.c ../svm.o svm_model_matlab.o

svmtrain312.$(MEX_EXT):       svmtrain312.c ../svm.h ../svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmtrain312.c ../svm.o svm_model_matlab.o

libsvmread312.$(MEX_EXT):	libsvmread312.c
	$(MEX) $(MEX_OPTION) libsvmread312.c

libsvmwrite312.$(MEX_EXT):	libsvmwrite312.c
	$(MEX) $(MEX_OPTION) libsvmwrite312.c

svm_model_matlab.o:     svm_model_matlab.c ../svm.h
	$(CXX) $(CFLAGS) -c svm_model_matlab.c

../svm.o: ../svm.cpp ../svm.h
	make -C .. svm.o

clean:
	rm -f *~ *.o *.mex* *.obj ../svm.o
