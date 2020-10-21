#include "mex.h"
#include "math.h"
#include "matrix.h"
//#include <omp.h>
#define SQRT_MAGIC_F 0x5f3759df
// 64  Bit float magic number
#define SQRT_MAGIC_D 0x5fe6ec85e7de30da

inline double fastSqrt_Bab(const double x)  
{
  union
  {
    long i;
    double x;
  } u;
  u.x = x;
  u.i = (((long)1)<<61) + (u.i >> 1) - (((long)1)<<51); 

  // One Babylonian Step
  u.x = 0.5F * (u.x + x/u.x);

  return u.x;
}

inline double invSqrt(const double x)
{
  const double xhalf = 0.5F*x;
 
  union // get bits for floating value
  {
    double x;
    long i;
  } u;
  u.x = x;
  u.i = SQRT_MAGIC_D - (u.i >> 1);  // gives initial guess y0
  return u.x*(1.5F - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
}

inline double fastSqrt_Q3(const double x)
{
  return x * invSqrt(x);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
int h, i, l, cnt, j;
// Define variables
int mrows1 	= mxGetM(prhs[0]);
int ncols1 	= mxGetN(prhs[0]);
int ngroups	= mxGetM(prhs[1]);
int ncols2	= mxGetN(prhs[1]);
int rowsB   	= ngroups+1;
int outputflag  = (int)mxGetScalar(prhs[2]);
int numthread   = (int)mxGetScalar(prhs[3]);
int chunk 	= ncols1/numthread;
double *L       = mxGetPr(prhs[1]);
double *Y;

if (nrhs == 5) { // Reorder matrix according to I (resampling)
	int mrowsI = mxGetM(prhs[4]);
	if (mrowsI > 0){
		//mexPrintf("\n\nLength of resampling vector = %d",mrows2);
		//mexPrintf("\n# of rows in Y = %d\n",mrows1);
		if (mrowsI != mrows1){
			mexErrMsgTxt("Length of resampling vector != rows of Y matrix!");	
		}
		double *oY = mxGetPr(prhs[0]);
		double *I = mxGetPr(prhs[4]); //Reorder index
		
		plhs[1] = mxCreateDoubleMatrix(mrows1,ncols1,mxREAL);
		Y = mxGetPr(plhs[1]);
		int ind;	
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
	else {
		Y = mxGetPr(prhs[0]);
	}
}
else if (nrhs == 4){ // Use input matrix as it is
	Y = mxGetPr(prhs[0]);
}
else { // Error!!!
	mexErrMsgTxt("Usage: pearson(<nxm data matrix>, <group boundary array>, <outputflag>, <no. of CPUs>, [resample index])");
}

if (ncols2 > 2)
mexErrMsgTxt("Label matrix has to have two columns and rows corresponding to the number of groups.");

int nbin = ngroups*(ngroups-1)/2;
int* bsamples = new int[rowsB];
int* nsamples = new int[ngroups];

// Output matrices
if (outputflag==0){
        plhs[0]=mxCreateDoubleMatrix(ncols1,1,mxREAL);
}
else if (outputflag == 1){
        plhs[0]=mxCreateDoubleMatrix(ncols1,nbin,mxREAL);
}
else {
	mexErrMsgTxt("Valid options for outputflag are [0=mean of binary comparisons],[1=binary comparisons].");
}
double *out=mxGetPr(plhs[0]);

bsamples[0] = 0;
for (h = 0; h<ngroups; ++h){
	bsamples[h+1] = L[ngroups + h] - 1;
	nsamples[h] = L[ngroups + h] - (L[h]-1);
}

int* nsamplesbin 	= new int[nbin];
double* avgLp 		= new double[nbin];
double* avgLn 		= new double[nbin];
double* dLv 		= new double[nbin];
double avgL, n1, n2;

cnt=0;
for (h = 0; h<ngroups-1; ++h)	
{
	n1=nsamples[h];
	for (l = h+1;l<ngroups;++l)
	{
		n2=nsamples[l];
		// Compute the group size for each binary comparison
		nsamplesbin[cnt] = n1+n2 ;
		
		// Now, compute the mean:
		avgL = (n1-n2) / (n1+n2);
		
		// ... and the positive and negative difference to the mean:
		avgLp[cnt] = 1-avgL;
		avgLn[cnt] = -1-avgL;
		
		// ... and the sum of square differences to the mean:
		dLv[cnt] = n1 * (avgLp[cnt] * avgLp[cnt]) + n2 * (avgLn[cnt] * avgLn[cnt]);
		++cnt;
	}
}	

/* Loop through all voxels */
#pragma omp parallel private(i,j,h,l,cnt) shared(Y,out,nsamplesbin,nsamples,bsamples,ngroups,rowsB,mrows1,ncols1,avgLp,avgLn,dLv,nbin)
{
double* sY =  new double[ngroups];
double pearsonsum, pearson, apearson;
double avgY, dY, dYv, s;
#pragma omp for schedule(static,chunk)
for (i = 0; i<ncols1; ++i)
{	
	/* Extract current LOO-data vector */
	//Compute sum of Y values for each group
	for (h=0;h<ngroups;++h)
	{
		sY[h]=0;
		for (j = bsamples[h]; j < bsamples[h+1]; ++j)
			sY[h] += Y[i*mrows1+j];
	}

	/* Loop through all groups */
	cnt=0; pearsonsum=0;
	for (h = 0; h<ngroups-1; ++h)	
	{
		for (l = h+1;l<ngroups;++l)
		{
			avgY = (sY[h]+sY[l]) / nsamplesbin[cnt];
			dYv = 0; s=0;
			for (j = bsamples[h]; j < bsamples[h+1]; ++j)
			{
				dY = Y[i*mrows1 + j] - avgY;
				dYv += dY*dY; /* Sum of square differences of Y */
				/* Sum of square differences of class vector [1 -1] has been computed above*/
				s += dY * avgLp[cnt]; 
			}
			for (j = bsamples[l]; j < bsamples[l+1]; ++j)
			{
				dY = Y[i*mrows1 + j] - avgY;
				dYv += dY*dY;
				s += dY * avgLn[cnt];
			}
			
			// Finally, compute correlation coefficient
			pearson = s/(sqrt(dYv)*sqrt(dLv[cnt]));
			apearson = fabs(pearson);
			pearsonsum += apearson;
			if (outputflag == 1){
				out[cnt*ncols1+i] = pearson;
			}
			++cnt;
		}
	}
	if (outputflag == 0){
		out[i] = pearsonsum/nbin;
    	}
}
delete [] sY;
}
delete [] bsamples;
delete [] nsamples;
delete [] nsamplesbin;
delete [] avgLp;
delete [] avgLn;
delete [] dLv;
}