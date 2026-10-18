/********************************************/
/*                                          */
/*Procedure de calcul de la densite uniforme*/
/*                                          */
/********************************************/

double calculLogDensiteUniformeSimplex(int nbIndividus,int nbDimensions,individu_t **individus);

/***********************************************/
/*                                             */
/*Procedure de maximisation de la vraisemblance*/
/*                                             */
/***********************************************/

int maximisationVraisemblanceSimplex(double **points,int nbDimensions,cluster_t cluster);

/************************************************************************************************/
/*                                                                                              */
/*Procedure de calcul de la log-vraisemblance d'une classe d'individus pour une loi de Dirichlet*/
/*                                                                                              */
/************************************************************************************************/

double calculLogVraisemblanceIndividusSimplex1cluster(int nbDimensions,double **points,
						      cluster_t cluster,double *alphas);

/********************************************************************/
/*                                                                  */
/*Procedure de calcul des log-vraisemblances pour tous les individus*/
/*par rapport a tous les clusters                                   */
/*                                                                  */
/********************************************************************/

void calculVraisemblancesPointsSimplex(int nbIndividus,int nbDimensions,individu_t **individus,
				       int nbClusters,cluster_t *clusters);
