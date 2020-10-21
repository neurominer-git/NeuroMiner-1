#include "mex.h"
#include "math.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double *out;
double *x2,*x1;
int i,j;
double avg1, avg2, var1, var2, avgdiff, val;
int mrows1,ncols1,mrows2,ncols2;

mrows1 	= mxGetM(prhs[0]);
ncols1 	= mxGetN(prhs[0]);
mrows2 	= mxGetM(prhs[1]);
ncols2 	= mxGetN(prhs[1]);
x1	= mxGetPr(prhs[0]);
x2	= mxGetPr(prhs[1]);

plhs[0]=mxCreateDoubleMatrix(ncols1,1,mxREAL);
out=mxGetPr(plhs[0]);

for (i = 0; i<ncols1; ++i)
{	
	avg1=0; avg2=0; var1=0; var2=0;
	
	/* First compute mean */
	for(j = 0; j<mrows1; ++j)
	{
		avg1= avg1 + x1[i*mrows1+j];
	}
	/* First compute mean */
	for(j = 0; j<mrows2; ++j)
	{
		avg2= avg2 + x2[i*mrows2+j];	
	}
	/* Check whether the current feature is empty */
	if (avg1 == 0 && avg2 == 0) continue;

	avg1 = avg1/mrows1;
	avg2 = avg2/mrows2;

	/* Then compute variance */
	for(j = 0; j<mrows1; ++j)
	{
		val = x1[i*mrows1+j] - avg1;
		var1 = var1 + val*val;
	}
	for(j = 0; j<mrows2; ++j)
	{
		val = x2[i*mrows2+j] - avg2;
		var2 = var2 + val*val;
	}
	var1 = var1/(mrows1-1);
	var2 = var2/(mrows2-1);

	/* Ok, now compute FDR for current feature */
	avgdiff = avg1 -avg2;
	out[i] = (avgdiff*avgdiff)/(var1+var2);
}
}
