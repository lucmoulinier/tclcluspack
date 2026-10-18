/**************************************************************/
/*                                                            */
/*Procedure de test pour savoir si l'on doit scinder un groupe*/
/*par la methode DPC (Density of Points Clustering)           */ 
/*                                                            */
/**************************************************************/

int testDPC(int nbIndividus,int nbDimensions,cluster_t *groupePere,
	    cluster_t *groupe1,cluster_t *groupe2,double **coordonneesIndividus,
	    double *cotesHyperpave,double *valeursMinimum,double risque,
	    int typeDonnees,int typeDensite,int donneesNormalisees,int nbSimulations);
