/********************************************************************/
/*                                                                  */
/*Procedure de recherche de coupe minimum normalisee de Xing et Karp*/
/*                                                                  */
/********************************************************************/

int normalizedCut(cluster_t *groupe,cluster_t *groupe1,cluster_t *groupe2);


/********************************************************/
/*                                                      */
/*Procedure renvoyant pour un vecteur donne la valeur   */
/*de la fonction objectif associee a la coupe normalisee*/
/*                                                      */
/********************************************************/

double evaluationCoupeNormaliseeMinimum(int nbDimensions,double *y,double *D,double **DminusW);
