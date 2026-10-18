/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/
/*le parametre diag indique si l'eventuelle matrice de variance covariance est diagonale*/
void kernelMixtureModels(int nbIndividus,
			 int nbDimensions,
			 individu_t **individus,
			 mixtureModel_t *mixtureModel,
			 int nbIterationsMax,
			 int nbClustersSelection,
			 int nbClustersMax,
			 int choixNoyau,
			 int choixDistribution,
			 int accelerationCalculs,
			 int nbPivots,
			 int methodeEstimationParametres,
			 double **kernels);

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*par modeles de melanges                                    */
/*                                                           */
/*************************************************************/

int kKernelMixtureModels(int nbIndividus,int nbDimensions,individu_t **individus,int nbClusters,
			 cluster_t *clusters,int nbIterationsMax,double *meilleureVraisemblance,
			 int choixNoyau,int choixDistribution,int accelerationCalculs,
			 int nbPivots,int methodeEstimationParametres,double **kernels);
