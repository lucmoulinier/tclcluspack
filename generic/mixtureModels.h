/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/

int clusteringMixtureModels(int nb_individus,
			    int nbDimensions,
			    individu_t **individus,
			    mixtureModel_t *mixturemodel,
			    int nb_iterations_max,
			    int nb_clusters_selection,
			    int type_donnees,
			    int diag,
			    int nbClustersMax,
			    int strictExecution);

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/

int kmixtureModels(int nb_individus,int nbDimensions,individu_t **individus,
		   int nb_clusters,cluster_t *clusters,int nb_iterations_max,
		   double *meilleure_vraisemblance,double logDensiteLoiUniforme,
		   int type_donnees,int diag,int strictExecution);

