#include "mex.h"
#include "math.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int h, i, j;
int mrows1 	= mxGetM(prhs[0]);
int ncols1 	= mxGetN(prhs[0]);
int ngroups	= mxGetM(prhs[1]);
int ncols2	= mxGetN(prhs[1]);
int numthread   = mxGetScalar(prhs[2]);
int chunk 	= ncols1/numthread;
double *L	= mxGetPr(prhs[1]);
double *Y;
if (nrhs == 4) { // Reorder matrix according to I (resampling)
	int mrowsI = mxGetM(prhs[3]);
		if (mrowsI > 0){
		//mexPrintf("\n\nLength of resampling vector = %d",mrows2);
		//mexPrintf("\n# of rows in Y = %d\n",mrows1);
		if (mrowsI != mrows1)
			mexErrMsgTxt("Length of resampling vector != rows of Y matrix!");	
		
		double *oY = mxGetPr(prhs[0]);
		double *I = mxGetPr(prhs[3]); //Reorder index
		int ind;
		plhs[1] = mxCreateDoubleMatrix(mrows1,ncols1,mxREAL);
		Y = mxGetPr(plhs[1]);
		#pragma omp parallel private(i,j, ind) shared(Y, I, oY, mrows1, ncols1)
		{
			#pragma omp for schedule(static,chunk)
			for (i = 0; i<ncols1; ++i){
				//mexPrintf("\nCol %d",i);
				for (j = 0;j<mrows1; ++j){
					ind = I[j]-1;
					//mexPrintf("\nRow %d => %d",j, ind+1);
					Y[i*mrows1+j] = oY[i*mrows1+ind];
				}
			}
		}
	}
}
else if (nrhs == 3){ // Use input matrix as it is
	Y = mxGetPr(prhs[0]);
}
else { // Error!!!
	mexErrMsgTxt("Usage: ftest(<nxm Y matrix>, <group boundary array>, <no. of CPUs>, [resample index]");
}

if (ncols2 > 2)
mexErrMsgTxt("Label matrix has to have two columns and rows corresponding to the number of groups.");

// Output matrix
plhs[0]=mxCreateDoubleMatrix(ncols1,1,mxREAL);
double *out=mxGetPr(plhs[0]);

double* nsamples =  new double[ngroups];
/* Compute group N */
for (h = 0; h<ngroups; ++h)	
	nsamples[h] = L[1*ngroups+h] - L[h];

#pragma omp parallel private(i,j,h) shared(Y, L, out, mrows1, nsamples, ngroups, ncols1)
{
double* avg =  new double[ngroups];
double* var =  new double[ngroups];
double val, grandavg, grandavgdiff, MSbetween, MSwithin;
/* Loop through all voxels */
#pragma omp for schedule(static,chunk)
for (i = 0; i<ncols1; ++i)
{	
	grandavg = 0; MSwithin=0; MSbetween=0;
	
	/* Loop through all groups */
	for (h = 0; h<ngroups; ++h)	
	{
		avg[h]=0; var[h]=0;
		
		/* First compute mean of current group*/
		for(j = L[h]-1; j<L[1*ngroups+h]-1; ++j)
		{
			avg[h] += Y[i*mrows1+j];
		}
		avg[h] /= nsamples[h];
		grandavg += avg[h];

		/* Now compute variance of current group */
		for(j = L[h]-1; j<L[1*ngroups+h]-1; ++j)
		{
			val = Y[i*mrows1+j] - avg[h];
			var[h] += val*val;
		}
		var[h] /= (nsamples[h]-1);
	}
	/* Finally compute the F-score*/
	/* F-score = Mean square of between variances / Mean square of within variances */
	grandavg /= ngroups;
	for (h = 0; h<ngroups; ++h)
	{
		grandavgdiff = avg[h]-grandavg;
		MSbetween += grandavgdiff*grandavgdiff;
		MSwithin += var[h];
	}
	MSbetween /= (ngroups-1);
	MSwithin /= ngroups;
	out[i] = MSbetween/MSwithin;
}
/* Free memory */
delete [] avg;
delete [] var;
}
delete [] nsamples;
}
