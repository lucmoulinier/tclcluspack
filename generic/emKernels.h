/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calculVraisemblancesKernels(int nbIndividus,int nbDimensions,individu_t **individus,int nbClusters,
				 cluster_t *clusters,int choixNoyau,int choixDistribution,
				 int methodeEstimationParametres,int accelerationCalculs,double **kernels);

/***********************************************************************/
/*                                                                     */
/*Procedure de calcul des parametres de la loi (exponentielle ou gamma)*/
/*qui renvoie la log-vraisemblance                                     */
/*                                                                     */
/***********************************************************************/

void estimationParametresKernelMoments(int nbDimensions,cluster_t *cluster,int accelerationCalculs,
				       int choixNoyau,int choixDistribution,double **kernels);

/*********************************************************/
/*                                                       */
/*Procedure d'ajustement de la constante de normalisation*/
/*                                                       */
/*********************************************************/

double normalisationKernel(int nbDimensions,int choixNoyau,int choixDistribution,cluster_t *cluster);

/***********************************************************************/
/*                                                                     */
/*Procedure de calcul des parametres de la loi (exponentielle ou gamma)*/
/*qui renvoie la log-vraisemblance                                     */
/*                                                                     */
/***********************************************************************/

void estimationParametresKernelMaximumLikelihood(int nbDimensions,cluster_t *cluster,int accelerationCalculs,
					         int choixNoyau,int choixDistribution,double **kernels,double regularisation);

/**************************************************/
/*                                                */
/*Calcul de la log vraisemblance pour la loi gamma*/
/*                                                */
/**************************************************/

double calculLogVraisemblanceGamma(int nbDimensions,cluster_t *cluster,int choixNoyau,double regularisation);
