/******************************************************************/
/*                                                                */
/*procedure d'apprentissage : k au sens des k plus proches voisins*/
/*                                                                */
/******************************************************************/

void apprentissageKernels(int *nbDimensions,int *tailleTrainingSet,individu_t **trainingSet,int tailleTestSet,
			  individu_t **testSet,int *nbClusters,cluster_t **clusters,mixtureModel_t ***mixtureModels,
			  mixtureModel_t *mixtureModelPourDensite,int bandwidthType,double bandwidth,int accelerationCalculs,
			  int *diag,int choixNoyau,int *choixDistribution,double scale,int methodeEstimationParametres,int verbose,
			  double regularisation);
