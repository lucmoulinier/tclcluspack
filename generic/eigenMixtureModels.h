/*************************************************/
/*                                               */
/*Clustering par eigen mixture models            */
/*Etape 1 : on surestime le nombre de composantes*/
/*Etape 2 : on fusionne les composantes par PCA  */
/*                                               */
/*************************************************/

void eigenMixtureModels(int nbIndividus,int nbDimensions,individu_t **individus,int nbIterationsMax,
			int nbClustersMax,int *nbClusters,int diag);

/***************************************************************/
/*                                                             */
/*Procedure de calcul de la pvaleur d'avoir la somme des carres*/
/*de k coordonnees d'un vecteur propre normalise superieure a x*/
/*                                                             */
/***************************************************************/

double pvaleursKcoordonneesVecteurPropre(int nbDimensions,int k,double x);

/*******************************************************/
/*                                                     */
/*Procedure de calcul de comparaison pour le classement*/
/*des valeurs propres de maniere croissante            */
/*                                                     */
/*******************************************************/

int compareCoordonnees(void const *a,void const *b);

/****************************************************************/
/*                                                              */
/*Procedure d'affectation des differentes composantes du melange*/
/*aux differents vecteurs propres                               */
/*                                                              */
/****************************************************************/

void affectationComposantesMelangeAuxVecteursPropres(int nbDimensions,double **vecteursPropres,
						     int *vecteursDesComposantes);
