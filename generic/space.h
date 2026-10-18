/**************************************************************/
/*                                                            */
/*Procedure de calcul du volume d'un tetraedre a nb_dimensions*/
/*                                                            */
/**************************************************************/

void calculVolumeSimplex(int nbDimensions,double *volumeSimplex);

/*************************************************************************************/
/*                                                                                   */
/*Procedure de calcul d'une sphere unitaire a nb_dimensions, dans R^{nb_dimensions+1}*/
/*                                                                                   */
/*************************************************************************************/

void calculSurfaceSphere(int nbDimensions,double *surfaceSphere);

/*********************************************************************************************/
/*                                                                                           */
/*Procedure de calcul de comparaison pour le classement des individus par distance croissante*/
/*                                                                                           */
/*********************************************************************************************/
int compareIndividus(void const *a,void const *b);

/********************************************************************************************/
/*                                                                                          */
/*Procedure de calcul de la densite non parametrique en utilisant les k plus proches voisins*/
/*                                                                                          */
/********************************************************************************************/

double calculDensiteNonParametrique(cluster_t *groupe,double *coordonnees,int dimension,
				    int nbVoisins,int ajout);

/*************************************************************/
/*                                                           */
/*Procedure de calcul des exposants d'Holder                 */
/*Methode : on forme un hyperpave avec pour longueur         */
/*la distances entre les deux coordonnees                    */
/*les plus eloignees sur chaque axe                          */
/*                                                           */
/*************************************************************/

void calculExposantsHolderPoints(int nbDimensions,int nbIndividusTotal,
				 double **coordonneesTousIndividus,int nbIndividus,
				 double **coordonneesIndividus,double *exposantsHolder,
				 double *cotesHyperpave,double *valeursMinimum,
				 int modeVerbose);

