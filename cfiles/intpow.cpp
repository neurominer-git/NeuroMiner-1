#include "mex.h"
#include "math.h"
#include "matrix.h"

double intpow( double base, int exponent ){
if (exponent == 0)
	return 1;
else
{
	double curpwr = base;
	int i;
	for(i=1; i<exponent; ++i)
		curpwr = base * curpwr;
	return curpwr;
}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double *out, pw;
double base = mxGetScalar(prhs[0]);
int exp = mxGetScalar(prhs[1]);
plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
out = mxGetPr(plhs[0]);
pw = intpow(base,exp);
mexPrintf("Power: %g", pw);
}