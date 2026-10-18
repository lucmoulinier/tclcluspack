/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/
/*le parametre diag indique si l'eventuelle matrice de variance covariance est diagonale*/

void GAPkmeans(int nbIndividus,int nbDimensions,individu_t **individus,int nbClustersMax,
	       int *meilleurNombreClusters);

/***************************************************************/
/*                                                             */
/*Procedure de calcul de l'inertie intra-classe d'un clustering*/
/*                                                             */
/***************************************************************/

double calculInertieIntraClasse(int nbIndividus,int nbDimensions,individu_t **individus,
				int nbClusters);
