# NeuroMiner-1
Prototype version of NeuroMiner used by the Section for Neurodiagnostic Applications (LMU) and developed by Prof. Nikolaos Koutsouleris

NeuroMiner has been developed to provide clinical researchers with cutting-edge machine learning methods for the analysis of heterogeneous data domains, such as clinical and neurocognitive read-outs, structural and functional neuroimaging data, and genetic information. The program can be considered as an interface to a large variety of unsupervised, and supervised pattern recognition algorithms that have been developed in the machine learning field over the last decades. Furthermore, the application implements itself different strategies for preprocessing, filtering and fusing heterogeneous data, training ensembles of predictors and visualizing and testing the significance of the computed predictive patterns. The current release candidate of NeuroMiner has been tested in the Section of Neurodiagnostic Applications on a variety of datasets from healthy controls and patients with psychiatric disorders and was designed specifically to create robust models with a high probability of generalization to new datasets. For reference, we include here a list of papers (see publications listed in https://www.ncbi.nlm.nih.gov/pubmed/?term=koutsouleris+n), which were all based on previous versions of the program.

More specifically, using a light-weight and interactive text-based menu system, NeuroMiner allows the user to:

* a) load their data easily (e.g., using spreadsheets, NifTi images, or SPM structures); 
* b) build a variety of cross-validation frameworks for classification and regression problems that have become a gold standard in the field (e.g., repeated nested cross-validation, leave-site-out cross-validation); 
* c) apply a range of preprocessing strategies (e.g., scaling, filtering, many forms of dimensionality reduction, etc.); 
* d) choose and combine cutting-edge supervised algorithms (e.g., support vector machine, elastic net, random forest, etc.); 
* e) apply feature selection procedures (e.g., wrappers), data fusion techniques, and stacked generalization;
* f) apply learned models to new data (external validation).  

To assist in selecting and analysing data, the user can visualise the data during input, monitor accuracy during learning, and understand the results of complex analyses using multiple display options. These allow the user to accurately report the data and also to understand the underlying machine learning analyses. Furthermore, the ability to apply the learned models to completely new data is important because it is quickly becoming a standard requirement of all machine learning studies. Combined, NeuroMiner gives the user the opportunity to design, implement, understand, and report machine learning analyses. 

DISCLAIMER

Please note that NeuroMiner is supplied as is and no formal maintenance is provided or implied. In no event shall the author of the software (heretofore known as the Author) be liable to any party for direct, indirect, special, incidental, or consequential damages, including lost profits, arising out of the use of this software and its documentation, even if the Author has been advised of the possibility of such damage. The Author specifically disclaims any warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. The software and accompanying documentation provided hereunder is provided 'as is'. The Author has no obligation to provide maintenance, support, updates, enhancements, or modifications (but we plan to). 
 
This is the beta release version of the software (NM 1; Elessar) and the software is undergoing regular updates. Please send any comments, questions, or bug reports to nm@pronia.eu or email.neurominer@gmail.com. 
