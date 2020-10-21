#ifndef _LIBSVM_H
#define _LIBSVM_H

#define LIBSVM_VERSION 289

#ifdef __cplusplus
extern "C" {
#endif

struct svm_performance{
	double tn;
	double tp;
	double fn;
	double fp;
	double fpr;
	double sens;
	double spec;
	double ppv;
	double npv;
	double mcc;
	double acc;
};
struct svm_decvalues{
	double n;
	double *decvalues;
	double *target;
};

/**
 * 	An element of a pattern in datasets
 * */
struct svm_node
{
	int index;		//	index of the item
	double value;	//	value of the item
};

/**
 * 	The problem: the dimension (l), the labels (*y) and the patterns (**x)
 * */
struct svm_problem
{
	int l;					//	number of elements in the input dataset
	double *y;				//	labels of the items in the input dataset
	struct svm_node **x;	//	items of the input dataset
};

/**
 * SVM Type: Classification, Classification with L2SVM, nu-Classification, 
 * One Class Classification, epsilon Regression, nu Regression, 
 * Support Vector Domain Description (SVDD), SVDD with L2SVM
 * 
 * C_SVC_L2, SVDD and SVDD_L2 were backported from the code available at
 * 
 *    http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/#18
 * 
 * written by Leland Wang and extended by Holger Froehlich.
 * 
 * --
 * Vincenzo Russo (vincenzo.russo@neminis.org)
 * */
enum { C_SVC, C_SVC_L2, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR, SVDD, SVDD_L2 };

/**
 * Kernel Type: Linear, Polynomial, Gaussian, Sigmoid, Stump, Perceptron, Laplace, Expo, Precomputed
 * 
 * Stump, Perceptron, Laplace, Expo were backported from the LIBSVM extension by Hsuan-Tien Lin
 * 
 * For his version of the software, see http://www.work.caltech.edu/~htlin/program/libsvm/#infensemble
 * For explanation about the four above kernels, see his Master's Thesis
 * 	"H.-T. Lin. Infinite Ensemble Learning with Support Vector Machines. 
 *   Master's Thesis, California Institute of Technology, May 2005."
 * http://www.work.caltech.edu/~htlin/publication/
 * 
 * --
 * Vincenzo Russo (vincenzo.russo@neminis.org)
 * */
enum { LINEAR, POLY, GAUSSIAN, SIGMOID, STUMP, PERC, LAPLACE, EXPO, PRECOMPUTED };	/* kernel_type */

/**
 * The parameters of a task
 * */
struct svm_parameter{
	int svm_type;		/* SVM Type */
	int kernel_type;	/* Kernel Type */
	int degree;			/* for poly */
	double gamma;		/* for poly/gaussian/sigmoid */
	double coef0;		/* for poly/sigmoid */

	/* these are for training only */
	double cache_size; 		/* in MB */
	double eps;				/* stopping criteria (the lower eps, the better the approximation) */
	double C;				/* for C_SVC, SVDD, EPSILON_SVR and NU_SVR */
	int nr_weight;			/* for C_SVC */
	int *weight_label;		/* for C_SVC */
	double* weight;			/* for C_SVC */
	double nu;				/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;				/* for EPSILON_SVR */
	int shrinking;			/* use the shrinking heuristics */
	int probability; 		/* do probability estimates */
	int decvals; /* decision values */
	int probability_fold; /* v-fold for probability value estimation */
};

/**
 * The model yielded by the training of an SVM
 * */
struct svm_model{
	struct svm_parameter param;	/* parameter */
	int nr_class;			/* number of classes, = 2 in regression/one class svm */
	int l;					/* total #SV*/
	struct svm_node **SV;			/* SVs (SV[l]) */
	
/* 	int *BSV_idx;			/* indices of the BSVs in the original dataset
 	int *SV_idx;			/* indices of the SVs in the original dataset
    int lbsv;				/* total #BSV */
	
	double **sv_coef;		/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */
	double *rho;			/* constants in decision functions (rho[k*(k-1)/2]) */
	double *w2;
	double *probA;          /* pariwise probability information */
	double *probB;
	/* for classification only */

	int *label;		/* label of each class (label[k])*/
	int *nSV;		/* number of SVs for each class (nSV[k])*/
					/* nSV[0] + nSV[1] + ... + nSV[k-1] = l*/
	int free_sv;		/* 1 if svm_model is created by svm_load_model*/
						/* 0 if svm_model is created by svm_train*/
};


/**
 * The decision function
 * */
struct decision_function{
	double *alpha;	/* Lagrangian multipliers */
	double rho;		
	double w2;
	
	int *BSV_idx;	/* indices of the BSVs in the original dataset*/
	int *SV_idx;	/* indices of the SVs in the original dataset*/
		
	int l,			/* number of total SVs (SVs + BSVs) */
		lbsv;		/* number of BSVs*/
};


struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param);
void svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, double *target);
void svm_cross_validation_dec(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, double *target, double *targetdec);

int svm_save_model(const char *model_file_name, const struct svm_model *model);
struct svm_model *svm_load_model(const char *model_file_name);

int svm_get_svm_type(const struct svm_model *model);
int svm_get_nr_class(const struct svm_model *model);
void svm_get_labels(const struct svm_model *model, int *label);
double svm_get_svr_probability(const struct svm_model *model);

void svm_predict_values(const struct svm_model *model, const struct svm_node *x, double* dec_values);
double svm_predict(const struct svm_model *model, const struct svm_node *x);
double svm_predict_probability(const struct svm_model *model, const struct svm_node *x, double* prob_estimates);

void svm_destroy_model(struct svm_model *model);
void svm_destroy_param(struct svm_parameter *param);

const char *svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);
int svm_check_probability_model(const struct svm_model *model);

extern void (*svm_print_string) (const char *);

#ifdef __cplusplus
}
#endif

#endif /* _LIBSVM_H */
