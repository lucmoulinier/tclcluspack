/******************************************************************/
/*                                                                */
/*procedure d'apprentissage : k au sens des k plus proches voisins*/
/*                                                                */
/******************************************************************/

void apprentissage(int *nbDimensions,int tailleTrainingSet,individu_t **trainingSet,int tailleTestSet,
		   individu_t **testSet,int *nbClusters,cluster_t **clusters,mixtureModel_t ***mixtureModel,
		   int *diag);
