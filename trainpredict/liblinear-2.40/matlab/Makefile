# This Makefile is used under Linux

MATLABDIR ?= /usr/local/matlab
CXX ?= g++
#CXX = g++-3.3
CC ?= gcc
CFLAGS = -Wall -Wconversion -O3 -fPIC -I$(MATLABDIR)/extern/include -I..

MEX = $(MATLABDIR)/bin/mex
MEX_OPTION = CC="$(CXX)" CXX="$(CXX)" CFLAGS="$(CFLAGS)" CXXFLAGS="$(CFLAGS)"
# comment the following line if you use MATLAB on a 32-bit computer
MEX_OPTION += -largeArrayDims
MEX_EXT = $(shell $(MATLABDIR)/bin/mexext)

all:	matlab

matlab:	binary

octave:
	@echo "please type make under Octave"
binary: train_liblin24.$(MEX_EXT) predict_liblin24.$(MEX_EXT) libsvmread_liblin24.$(MEX_EXT) libsvmwrite_liblin24.$(MEX_EXT)

train_liblin24.$(MEX_EXT): train_liblin24.c ../linear.h ../newton.cpp ../linear.cpp linear_model_matlab.c \
		../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c
	$(MEX) $(MEX_OPTION) train_liblin24.c ../newton.cpp ../linear.cpp linear_model_matlab.c \
		../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c

predict_liblin24.$(MEX_EXT): predict_liblin24.c ../linear.h ../newton.cpp ../linear.cpp linear_model_matlab.c \
		../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c
	$(MEX) $(MEX_OPTION) predict_liblin24.c ../newton.cpp ../linear.cpp linear_model_matlab.c \
		../blas/daxpy.c ../blas/ddot.c ../blas/dnrm2.c ../blas/dscal.c

libsvmread_liblin24.$(MEX_EXT):	libsvmread_liblin24.c
	$(MEX) $(MEX_OPTION) libsvmread_liblin24.c

libsvmwrite_liblin24.$(MEX_EXT):	libsvmwrite_liblin24.c
	$(MEX) $(MEX_OPTION) libsvmwrite_liblin24.c

clean:
	rm -f *~ *.o *.mex* *.obj
