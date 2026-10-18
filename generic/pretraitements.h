
/*************************************/
/*                                   */
/*Procedure de reduction de dimension*/
/*                                   */
/*************************************/

void reductionDimension(int nbIndividus,individu_t **individus,int nbDimensions,
			int nbDimensionsDemandees);

/********************************************/
/*                                          */
/*Procedure de standardisation des individus*/
/*                                          */
/********************************************/

/* renvoie -1 si tout va bien, sinon l'indice de l'individu dont l'ecart-type est nul */
int standardisationIndividus(int nbIndividus,int nbDimensions,individu_t **individus);

/**************************************************/
/*                                                */
/*Procedure de projection des points sur la sphere*/
/*                                                */
/**************************************************/

void projectionSurSphere(int nbIndividus,int nbDimensions,individu_t **individus);

/**************************************************************/
/*                                                            */
/*Procedure de calcul des similarites a partir des coordonnees*/
/*en utilisant les correlations                               */
/*                                                            */
/**************************************************************/

void coordonneesToSimilarites(int nbIndividus,int nbDimensions,individu_t **individus);

/*********************************************************************/
/*                                                                   */
/*Procedure de reduction en utilisant un noyau                       */
/*on utilise des individus parmi individus1 pour reduire la dimension*/
/*                                                                   */
/*********************************************************************/

void reductionDimensionParNoyau(int nbIndividus1,int nbDimensions,individu_t **individus1,int nbIndividus2,
				individu_t **individus2,int choixNoyau,int nbPivots,int *newNbDimensions);
