#define GROUPE1 0
#define GROUPE2 1

/******************************************************/
/*                                                    */
/*Procedure de realisation d'un clustering d'un groupe*/ 
/*en deux par la methode des k-means avec k=2         */
/*                                                    */
/******************************************************/

int twoMeans(int nbDimensions,cluster_t *groupe,cluster_t *groupe1,cluster_t *groupe2,
	      int typeDonnees);

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectationIndividus2means(int nbDimensions,int nbIndividus,individu_t **individus,
				double *centre1,double *centre2,int *groupesDesIndividus,
				double *inertieIntraClasse,int typeDonnees,double *centre);
