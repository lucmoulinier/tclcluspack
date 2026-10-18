/********************************************/
/*                                          */
/*Procedure de calcul de la densite uniforme*/
/*                                          */
/********************************************/

double calculLogDensiteUniformeGaussienneSpherique(int nbIndividus,int nbDimensions,individu_t **individus);

/***********************************************/
/*                                             */
/*Procedure de maximisation de la vraisemblance*/
/*                                             */
/***********************************************/

void maximisationVraisemblanceGaussienneSpherique(int nbPoints,double **points,int nb_dimensions,
						  cluster_t cluster);

/**************************************************/
/*                                                */
/*Procedure de calcul du mu moyen pour initialiser*/
/*la maximisation de la vraisemblance             */
/*                                                */
/**************************************************/

void calculMuMoyenGaussienneSpherique(int nbPoints,double **coordonnees,int nbDimensions,double *mu);

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calculVraisemblancesPointsGaussienneSpherique(int nbIndividus,int nbDimensions,individu_t **individus,
						   int nbClusters,cluster_t *clusters);

/***********************************************/
/*                                             */
/*Procedure de calcul des nouvelles coordonnees*/
/*apres rotation du repere                     */
/*                                             */
/***********************************************/

void calculCoordonneesApresRotation(int nbDimensions,int nbPoints,double **points,
				    double *angles,double **pointsRot);

/*********************************************************/
/*                                                       */
/*Procedure de calcul des coordonnees des points projetes*/
/*                                                       */
/*********************************************************/

void calculCoordonneesApresProjection(int nbDimensions,int nbPoints,double **pointsRot,
				      double **pointsProjetes);

/**********************************************/
/*                                            */
/*Procedure de calcul des coordonnees polaires*/
/*en fonction des coordonnees cartesiennes    */
/*                                            */
/**********************************************/

void cartesiennes2polaires(int nbDimensions,double *cart,double *pol);
