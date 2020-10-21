#include "mex.h"
#include "math.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double *oY 	= mxGetPr(prhs[0]);
double *I 	= mxGetPr(prhs[1]); //Reorder index
int mrows1 	= mxGetM(prhs[0]);
int ncols1 	= mxGetN(prhs[0]);
int mrowsI 	= mxGetM(prhs[1]);
int numthreads = mxGetScalar(prhs[2]);
int chunk = ncols1/numthreads;
int i, j, ind;
if (mrowsI != mrows1)
	mexErrMsgTxt("Length of resampling vector ~= rows of Y matrix!");	

plhs[0] = mxCreateDoubleMatrix(mrows1,ncols1,mxREAL);
double *Y = mxGetPr(plhs[0]);

#pragma omp parallel private(i,j,ind) shared(Y, I, oY, mrows1, ncols1)
{
	#pragma omp for schedule(static,chunk)
	for (i = 0; i<ncols1; ++i){
		//mexPrintf("\nCol %d",i);
		for (j = 0;j<mrows1; ++j){
			ind = (int)I[j]-1;
			//mexPrintf("\nRow %d => %d",j, ind+1);
			Y[i*mrows1+j] = oY[i*mrows1+ind];
		}
	}
}

}