/********************************************/
/*                                          */
/*Procedure de calcul de la densite uniforme*/
/*                                          */
/********************************************/

double calculLogDensiteUniformeGaussienne(int nbIndividus,int nbDimensions,individu_t **individus);

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calculVraisemblancesPointsGaussienne(int nbIndividus,int nbDimensions,individu_t **individus,
					     int nbClusters,cluster_t *clusters,int diag);

/**************************************************/
/*                                                */
/*Procedure de calcul du mu moyen pour initialiser*/
/*la maximisation de la vraisemblance             */
/*                                                */
/**************************************************/

void calculMuMoyen(int nbPoints,double **coordonnees,int nbDimensions,double *mu);

/****************************************************************************/
/*                                                                          */
/*Procedure de maximisation de la vraisemblance pour des donnees gaussiennes*/
/*                                                                          */
/****************************************************************************/
  
int maximisationVraisemblanceGaussienne(double **points,int nbDimensions,cluster_t cluster,
					  int diag);
