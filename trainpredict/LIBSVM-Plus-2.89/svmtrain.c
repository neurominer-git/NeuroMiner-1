#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "svm.h"
#include <math.h>
#include "mex.h"
#include "svm_model_matlab.h"
#include "matrix.h"

#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif 

#define CMD_LEN 2048
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))

void print_null(const char *s) {}

void exit_with_help()
{
	mexPrintf(
	"Usage: model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
	"options:\n"
	"-s svm_type : set type of SVM (default 0)\n"
	"\t0 -- C-SVC (L1SVM)\n"
	"\t1 -- C-SVC (L2SVM)\n"
	"\t2 -- nu-SVC\n"
	"\t3 -- one-class SVM\n"
	"\t4 -- epsilon-SVR\n"
	"\t5 -- nu-SVR\n"
	"\t6 -- SVDD (L1SVM)\n"
	"\t7 -- SVDD (L2SVM)\n"
	"-t kernel_type : set type of kernel function (default 2)\n"
	"\t0 -- linear: u'*v\n"
	"\t1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	"\t2 -- gaussian: exp(-gamma*||u-v||^2)\n"
	"\t3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	"\t4 -- stump: -|u-v| + coef0\n"
	"\t5 -- perceptron: -||u-v|| + coef0\n"
	"\t6 -- laplacian: exp(-gamma*|u-v|)\n"
	"\t7 -- exponential: exp(-gamma*||u-v||)\n"
	"\t8 -- precomputed kernel (kernel values in training_set_file)\n"
	"-d degree: set degree in kernel function (default 3)\n"
	"-g gamma: set gamma in kernel function (default 1/k)\n"
	"-r coef0: set coef0 in kernel function (default 0)\n"
	"-c cost: set the parameter C of C-SVC (L1/L2-SVM), SVDD (L1/L2-SVM), epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu: set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon: set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-m cachesize: set cache memory size in MB (default 100)\n"
	"-e epsilon: set tolerance of termination criterion (default 0.001)\n"
	"-h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	"-b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
	"-wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
	"-v n: n-fold cross validation mode\n"
	"-q : quiet mode (no outputs)\n"
	"-u n: n-fold for probability value estimation (works only with -b 1)"
	);
}
/*
void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}*/

// void parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name);
// void read_problem(const char *filename);
// void do_cross_validation();

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
// struct svm_performance svmperf;
// struct svm_decvalues svmdec;
// struct svm_performance mPerf;
// struct svm_decvalues mDec;
//////
void do_cross_validation(struct svm_performance*, struct svm_decvalues*);

/**
 * storeValue( ) - store the given double value in the mxStruct
 * structure at field fieldName */
void storeValue(mxArray *mxStruct, char *fieldName, double value)
{
  char msg[1025];
  mxArray *fieldValue; 

  int fieldNum = mxGetFieldNumber(mxStruct, fieldName); 
  if (fieldNum < 0) {
    sprintf(msg,"Error: unknown field: (%s)\n", fieldName);
    mexErrMsgTxt(msg);
  }

  fieldValue = mxCreateDoubleMatrix(1,1,mxREAL); 
  *mxGetPr(fieldValue) = value; 
  mxSetFieldByNumber(mxStruct, 0, fieldNum, fieldValue); 
//mxDestroyArray(fieldValue);
}

/**
 * storeArray( ) - store the given array double array with n elements
 * into the structure mxStruct at field named name.
 */
void storeArray(mxArray *mxStruct, char *name, double *srcArray, int n)
{
  mxArray *fieldValue;
  double *tmpArray;
  int i;

  fieldValue = mxCreateDoubleMatrix(n,1,mxREAL); 
  tmpArray = mxGetPr(fieldValue);

  for (i = 0; i < n; i++) 
    tmpArray[i] = srcArray[i];

  mxSetField(mxStruct, 0, name, fieldValue); 
//mxDestroyArray(fieldValue);
}

// void do_cross_validation(struct svm_performance);

int cross_validation;
int nr_fold;

void do_cross_validation(mxArray *mPtr)
{
	int i;
	int total_correct = 0;
	double tp = 0;
	double tn = 0;
	double fp = 0;
	double fn = 0;
// 	double sens = 0.0;
// 	double spec = 0.0;
// 	double fpr = 0.0;
// 	double ppv = 0.0;
// 	double npv = 0.0;
// 	double mcc = 0.0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double,prob.l);
	double retval = 0.0;
	double *targetdec = Malloc(double,prob.l);
	if (param.decvals == 1)
	{
		svm_cross_validation_dec(&prob,&param,nr_fold,target,targetdec);
		storeArray(mPtr, "dec_values", targetdec, prob.l);
		storeArray(mPtr, "target", target, prob.l);
// 		for(i=0;i<prob.l;i++)
// 		{
// 			decval[i] = targetdec[i];
// 			mPtr2->target[i] = target[i];
// 		}*/

// 		mPtr2->n = (double)prob.l;
	}
	else
	{	
		svm_cross_validation(&prob,&param,nr_fold,target);
		storeValue(mPtr, "dec_values", 0);
		storeValue(mPtr, "target", 0);
	}

	if(param.svm_type == EPSILON_SVR ||
	   param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		mexPrintf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		mexPrintf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
		retval = total_error/prob.l;
	}
	else
	{
		for(i=0;i<prob.l;i++)
			{
// 			mexPrintf("Target[i] = %2.2f, prob.y[i] = %g\n",target[i],prob.y[i]);
			if(target[i] == prob.y[i])
			{
				++total_correct;
			}
			if (target[i] < 0 && prob.y[i] < 0)
			{
				++tn;
			}
			else if (target[i] > 0 && prob.y[i] < 0)
			{
				++fp;
			}
			else if (target[i] > 0 && prob.y[i] > 0)
			{
				++tp;
			}
			else if (target[i] < 0 && prob.y[i] > 0)
			{
				++fn;
			}
		}
		storeValue(mPtr, "TP", tp);
		storeValue(mPtr, "FP", fp);
		storeValue(mPtr, "TN", tn);
		storeValue(mPtr, "FN", tp);
		storeValue(mPtr, "FPR", fp / (fp + tn)*100.0);
		storeValue(mPtr, "sens", tp / (tp + fn)*100.0);
		storeValue(mPtr, "spec", tn / (tn + fp)*100.0);
		storeValue(mPtr, "PPV", tp / (tp + fp)*100.0);
		storeValue(mPtr, "NPV", tn / (tn + fn)*100.0);
		storeValue(mPtr, "MCC", (tp*tn-fp*fn)/sqrt((tp+fn)*(fp+tn)*(tp+fp)*(fn+tn)));
		storeValue(mPtr, "acc", 100.0*total_correct/prob.l);
// 		mexPrintf("TP, TN, FP, FN = %2.0f, %2.0f, %2.0f, %2.0f\n",mPtr->tp,mPtr->tn,mPtr->fp,mPtr->fn);
// 		mexPrintf("Cross Validation Accuracy = %g%%\n",mPtr->acc);
// 		mexPrintf("Sensitivity = %g%%\n",mPtr->sens);
// 		mexPrintf("Specificity = %g%%\n",mPtr->spec);
// 		mexPrintf("PPV = %g%%\n",mPtr->ppv);
// 		mexPrintf("NPV = %g%%\n",mPtr->npv);
// 		mexPrintf("MCC = %1.4f\n",mPtr->mcc);
	}
	free(target);
	free(targetdec);
}

// nrhs should be 3
int parse_command_line(int nrhs, const mxArray *prhs[], char *model_file_name)
{
	int i, argc = 1;
	char cmd[CMD_LEN];
	char *argv[CMD_LEN/2];

	// default values
	param.svm_type = C_SVC;
	param.kernel_type = GAUSSIAN;
	param.degree = 3;
	param.gamma = 0;	// 1/k
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	param.decvals = 0;
	param.probability_fold = 5;
	cross_validation = 0;

	if(nrhs <= 1)
		return 1;

	if(nrhs > 2)
	{
		// put options in argv[]
		mxGetString(prhs[2], cmd, mxGetN(prhs[2]) + 1);
		if((argv[argc] = strtok(cmd, " ")) != NULL)
			while((argv[++argc] = strtok(NULL, " ")) != NULL)
				;
	}

	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			return 1;
		switch(argv[i-1][1])
		{
			case 's':
				param.svm_type = atoi(argv[i]);
				break;
			case 't':
				param.kernel_type = atoi(argv[i]);
				break;
			case 'd':
				param.degree = atoi(argv[i]);
				break;
			case 'g':
				param.gamma = atof(argv[i]);
				break;
			case 'r':
				param.coef0 = atof(argv[i]);
				break;
			case 'n':
				param.nu = atof(argv[i]);
				break;
			case 'm':
				param.cache_size = atof(argv[i]);
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'e':
				param.eps = atof(argv[i]);
				break;
			case 'p':
				param.p = atof(argv[i]);
				break;
			case 'h':
				param.shrinking = atoi(argv[i]);
				break;
			case 'b':
				param.probability = atoi(argv[i]);
				break;
			case 'q':
				svm_print_string = &print_null;
				i--;
				break;
			case 'v':
				cross_validation = 1;
				nr_fold = atoi(argv[i]);
				if(nr_fold < 2)
				{
					mexPrintf("n-fold cross validation: n must >= 2\n");
					return 1;
				}
				break;
			case 'w':
				++param.nr_weight;
				param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
				param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
				param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
				param.weight[param.nr_weight-1] = atof(argv[i]);
				break;
			case 'a':
				param.decvals = atoi(argv[i]);
				break;
			case 'u':
				param.probability_fold = atoi(argv[i]);
				break;
			default:
				mexPrintf("unknown option\n");
				return 1;

		}
	}
	return 0;
}

// read in a problem (in svmlight format)
int read_problem_dense(const mxArray *label_vec, const mxArray *instance_mat)
{
	int i, j, k;
	int elements, max_index, sc, label_vector_row_num;
	double *samples, *labels;

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;
	
	labels = mxGetPr(label_vec);
	samples = mxGetPr(instance_mat);
	sc = mxGetN(instance_mat);

	elements = 0;
	// the number of instance
	prob.l = mxGetM(instance_mat);
	label_vector_row_num = mxGetM(label_vec);

	if(label_vector_row_num!=prob.l)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		return -1;
	}

	if(param.kernel_type == PRECOMPUTED)
		elements = prob.l * (sc + 1);
	else
	{
		for(i = 0; i < prob.l; i++)
		{
			for(k = 0; k < sc; k++)
				if(samples[k * prob.l + i] != 0)
					elements++;
			// count the '-1' element
			elements++;
		}
	}

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	max_index = sc;
	j = 0;
	for(i = 0; i < prob.l; i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];

		for(k = 0; k < sc; k++)
		{
			if(param.kernel_type == PRECOMPUTED || samples[k * prob.l + i] != 0)
			{
				x_space[j].index = k + 1;
				x_space[j].value = samples[k * prob.l + i];
				j++;
			}
		}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				mexPrintf("Wrong input format: sample_serial_number out of range\n");
				return -1;
			}
		}

	return 0;
}

int read_problem_sparse(const mxArray *label_vec, const mxArray *instance_mat)
{
	int i, j, k, low, high;
	mwIndex *ir, *jc;
	int elements, max_index, num_samples, label_vector_row_num;
	double *samples, *labels;
	mxArray *instance_mat_col; // transposed instance sparse matrix

	prob.x = NULL;
	prob.y = NULL;
	x_space = NULL;

	// transpose instance matrix
	{
		mxArray *prhs[1], *plhs[1];
		prhs[0] = mxDuplicateArray(instance_mat);
		if(mexCallMATLAB(1, plhs, 1, prhs, "transpose"))
		{
			mexPrintf("Error: cannot transpose training instance matrix\n");
			return -1;
		}
		instance_mat_col = plhs[0];
		mxDestroyArray(prhs[0]);
	}

	// each column is one instance
	labels = mxGetPr(label_vec);
	samples = mxGetPr(instance_mat_col);
	ir = mxGetIr(instance_mat_col);
	jc = mxGetJc(instance_mat_col);

	num_samples = mxGetNzmax(instance_mat_col);

	// the number of instance
	prob.l = mxGetN(instance_mat_col);
	label_vector_row_num = mxGetM(label_vec);

	if(label_vector_row_num!=prob.l)
	{
		mexPrintf("Length of label vector does not match # of instances.\n");
		return -1;
	}

	elements = num_samples + prob.l;
	max_index = mxGetM(instance_mat_col);

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node, elements);

	j = 0;
	for(i=0;i<prob.l;i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];
		low = jc[i], high = jc[i+1];
		for(k=low;k<high;k++)
		{
			x_space[j].index = ir[k] + 1;
			x_space[j].value = samples[k];
			j++;
	 	}
		x_space[j++].index = -1;
	}

	if(param.gamma == 0)
		param.gamma = 1.0/max_index;

	return 0;
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

// Interface function of matlab
// now assume prhs[0]: label prhs[1]: features
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
{
	const char *error_msg;
// 	int i;

	// fix random seed to have same results for each run
	// (for cross validation and probability estimation)
	srand(1);

	// Transform the input Matrix to libsvm format
	if(nrhs > 0 && nrhs < 4)
	{
		int err;


		if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
			mexPrintf("Error: label vector and instance matrix must be double\n");
			fake_answer(plhs);
			return;
		}

		if(parse_command_line(nrhs, prhs, NULL))
		{
			exit_with_help();
			svm_destroy_param(&param);
			fake_answer(plhs);
			return;
		}

		if(mxIsSparse(prhs[1]))
		{
			if(param.kernel_type == PRECOMPUTED)
			{
				// precomputed kernel requires dense matrix, so we make one
				mxArray *rhs[1], *lhs[1];

				rhs[0] = mxDuplicateArray(prhs[1]);
				if(mexCallMATLAB(1, lhs, 1, rhs, "full"))
				{
					mexPrintf("Error: cannot generate a full training instance matrix\n");
					svm_destroy_param(&param);
					fake_answer(plhs);
					return;
				}
				err = read_problem_dense(prhs[0], lhs[0]);
				mxDestroyArray(lhs[0]);
				mxDestroyArray(rhs[0]);
			}
			else
				err = read_problem_sparse(prhs[0], prhs[1]);
		}
		else
			err = read_problem_dense(prhs[0], prhs[1]);

		// svmtrain's original code
		error_msg = svm_check_parameter(&prob, &param);

		if(err || error_msg)
		{
			if (error_msg != NULL)
				mexPrintf("Error: %s\n", error_msg);
			svm_destroy_param(&param);
			free(prob.y);
			free(prob.x);
			free(x_space);
			fake_answer(plhs);
			return;
		}

		if(cross_validation)
		{
// 			double *ptr;
			const char *field_names[] = { "TP", "FP", "TN", "FN", "sens", "spec", "FPR", "PPV", "NPV", "MCC", "acc","dec_values","target" };
			mxArray *mxStruct;
  			mxStruct = mxCreateStructMatrix(1,1,NUMBER_OF_FIELDS,field_names);
			
			do_cross_validation(mxStruct);

// 			storeValue(mxStruct, "TP", mf->tp);
//   			storeValue(mxStruct, "FP", mf->fp);
//   			storeValue(mxStruct, "TN", mf->tn);
//   			storeValue(mxStruct, "FN", mf->fn);
//   			storeValue(mxStruct, "sens", mf->sens);
// 			storeValue(mxStruct, "spec", mf->spec);
// 			storeValue(mxStruct, "FPR", mf->fpr);
// 			storeValue(mxStruct, "PPV", mf->ppv);
// 			storeValue(mxStruct, "NPV", mf->npv);
// 			storeValue(mxStruct, "MCC", mf->mcc);
// 			storeValue(mxStruct, "acc", mf->acc);

// 			ptr = mxGetPr(plhs[0]);
			plhs[0] = mxStruct;
			//mxDestroyArray(mxStruct);
			
		}
		else
		{
			int nr_feat = mxGetN(prhs[1]);
			const char *error_msg;
			model = svm_train(&prob, &param);
			error_msg = model_to_matlab_structure(plhs, nr_feat, model);
			if(error_msg)
				mexPrintf("Error: can't convert libsvm model to matrix structure: %s\n", error_msg);
			svm_destroy_model(model);
		}
		svm_destroy_param(&param);
		free(prob.y);
		free(prob.x);
		free(x_space);
		return;
	}
	else
	{
		exit_with_help();
		fake_answer(plhs);
		return;
	}
}

