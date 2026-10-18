/*****************************************/
/*                                       */
/*Procedure de prediction sur le test set*/
/*                                       */
/*****************************************/

void predictionKernels(int nbDimensions,int nbClusters,cluster_t *clusters,int tailleTrainingSet,individu_t **trainingSet,
		       int tailleTestSet,individu_t **testSet,
			mixtureModel_t **mixtureModels,mixtureModel_t mixtureModelPourDensite,int bandwidthType,
		       double bandwidth,int accelerationCalculs,int diag,int choixNoyau,int choixDistribution,
		       int methodeEstimationParametres);

