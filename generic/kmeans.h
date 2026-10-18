/****************************************************/
/*                                                  */
/*Procedure de realisation de clustering par k-means*/
/*en partant d'un debut de clustering               */ 
/*                                                  */
/****************************************************/

int kmeans(int nbIndividus,int nbDimensions,individu_t **individus,int nbClasses);

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectationIndividusKmeans(int nbDimensions,int nbIndividus,int nbClasses,
			       int *nbIndividusParClasse,individu_t **individus,
			       double **centres,double *inertieIntraClasse,double *centreGravite);


