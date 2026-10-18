
/********************************************************************************************/
/*                                                                                          */
/*Procedure de recherche des k plus proches voisins pour un cluster avec stockage des pivots*/
/*                                                                                          */
/********************************************************************************************/

void kNearestNeighbours(int nbPivots,int kPlusProches,int nbIndividus,int nbDimensions,individu_t **individus,
			cluster_t *cluster);

/*****************************************************************************************/
/*                                                                                       */
/*Procedure de recherche des k plus proches voisins pour un individu dans un cluster*/
/*                                                                                       */
/*****************************************************************************************/

void kNearestNeighboursInOtherCluster(int nbPivots,int kPlusProches,int nbDimensions,individu_t *individu,
				      cluster_t *cluster);
