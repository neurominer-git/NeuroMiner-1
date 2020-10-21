# This Makefile is used under Linux

MATLABDIR = /lrz/sys/applications/matlab/R2011a/
#MATLABDIR = /sbia/sbiasfw/external/matlab/R2009B/
#MATLABDIR = /MRprogs3/MATLAB/R2011a/
# for Mac
# MATLABDIR ?= /opt/local/matlab

CXX ?= g++
#CXX = g++-4.1
CFLAGS = -Wall -Wconversion -O3 -fPIC -I$(MATLABDIR)/extern/include

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC\#$(CXX) CXX\#$(CXX) CFLAGS\#"$(CFLAGS)" CXXFLAGS\#"$(CFLAGS)"
# comment the following line if you use MATLAB on 32-bit computer
MEX_OPTION += -largeArrayDims
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)

OCTAVEDIR ?= /usr/include/octave
OCTAVE_MEX = env CC=$(CXX) mkoctfile
OCTAVE_MEX_OPTION = --mex
OCTAVE_MEX_EXT = mex
OCTAVE_CFLAGS = -Wall -O3 -fPIC -I$(OCTAVEDIR)

all:	matlab

matlab:	binary

octave:
	@make MEX="$(OCTAVE_MEX)" MEX_OPTION="$(OCTAVE_MEX_OPTION)" \
	MEX_EXT="$(OCTAVE_MEX_EXT)" CFLAGS="$(OCTAVE_CFLAGS)" \
	binary

binary: svmpredict291.$(MEX_EXT) svmtrain291.$(MEX_EXT) libsvmread291.$(MEX_EXT) libsvmwrite291.$(MEX_EXT)

svmpredict291.$(MEX_EXT):     svmpredict291.c svm.h svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmpredict291.c svm.o svm_model_matlab.o

svmtrain291.$(MEX_EXT):       svmtrain291.c svm.h svm.o svm_model_matlab.o
	$(MEX) $(MEX_OPTION) svmtrain291.c svm.o svm_model_matlab.o

libsvmread291.$(MEX_EXT):	libsvmread291.c
	$(MEX) $(MEX_OPTION) libsvmread291.c

libsvmwrite291.$(MEX_EXT):	libsvmwrite291.c
	$(MEX) $(MEX_OPTION) libsvmwrite291.c

svm_model_matlab.o:     svm_model_matlab.c svm.h
	$(CXX) $(CFLAGS) -c svm_model_matlab.c

svm.o:  svm.cpp svm.h
	$(CXX) $(CFLAGS) -c svm.cpp

clean:
	rm -f *~ *.o *.mex* *.obj
