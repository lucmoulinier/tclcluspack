/****************************************************/
/*                                                  */
/*Procedure de realisation de clustering par k-means*/
/*en partant d'un debut de clustering               */ 
/*                                                  */
/****************************************************/
/*typeNoyau peut valoir AUCUN_NOYAU ou RADIAL_GAUSSIAN*/
int kernelKmeans(int nbIndividus,
		 int nbDimensions,
		 individu_t **individus,
		 int nbClusters,
		 int choixNoyau,
		 int accelerationCalculs,
		 double **kernels);

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*pour la methode de kernel kmeans                 */
/*                                                 */
/***************************************************/

int affectationIndividusKernelKmeans(int nbDimensions,int nbIndividus,int nbClusters,
				     cluster_t *clusters,individu_t **individus,
				     double *inertieIntraClasses,int initialisation,
				     int choixNoyau,int accelerationCalculs,double **kernels);
