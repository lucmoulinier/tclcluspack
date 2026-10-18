/******************************************************/
/*                                                    */
/*Procedure de filtrage des dimensions ou celles ayant*/ 
/*des variances sont trop petites sont eliminees      */
/*                                                    */
/******************************************************/

void filtreDimensions(int nbIndividus,int nbDimensionsOrigine,int *nbDimensionsUtilisees,
		      individu_t **individus);


/******************************************************/
/*                                                    */
/*Procedure de filtrage des dimensions ou celles ayant*/ 
/*des variances sont trop petites sont eliminees      */
/*                                                    */
/******************************************************/

void filtreDimensionsDeuxEnsembles(int nbIndividus,int nbDimensionsOrigine,int *nbDimensionsUtilisees,
				   individu_t **individus,int nbIndividus2,individu_t **individus2);

/*********************************************/
/*                                           */
/*Procedure de filtrage des individus suivant*/ 
/*une valeur seuil et le type des donnees    */
/*                                           */
/*********************************************/

void filtreIndividus(int *nbIndividus,int nbDimensions,individu_t **individus,
		     int typeDonnees,double distanceMin,int distancesDirectes);
