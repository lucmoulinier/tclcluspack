/********************************************************************************/
/*                                                                              */
/*Procedure de decoupage d'arbre avec le test anova sur les densites de presence*/
/*                                                                              */
/********************************************************************************/

void decoupageArbreDPC(int nbIndividus,int nbDimensions,individu_t **individus,
		       noeud_t *racine,int *nbClusters,int typeDonnees,
		       int nbClustersSelection,int donneesNormalisees,int typeDensite,
		       int nbVoisins);
