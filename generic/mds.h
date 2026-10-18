/***************************************************************************************/
/*                                                                                     */
/*Procedure de calcul des matrices A et B obtenues a partir de la matrice des distances*/
/*                                                                                     */
/***************************************************************************************/     

void matriceAB(int nbIndividus,double **individus,double **A,double **B);


/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des coordonnees normalisees qui viennent de MDS*/
/*                                                                   */
/*********************************************************************/

void distancesToCoordonnees(int nbIndividus,individu_t **individus,double **distances,
			    int nbDimensions);
