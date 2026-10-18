/*choix distribution*/
#define EXPONENTIAL 1
#define GAMMA 2
#define BETA 3

/*constantes pour les noyaux*/
#define AUCUN_NOYAU 0
#define RADIAL_GAUSSIAN 1
#define EPANECHNIKOV 2
#define POLYNOMIAL2 3
#define POLYNOMIAL3 4
#define POLYNOMIAL4 5

/*constantes pour les bandwidths*/
#define FIXED_BANDWIDTH 0
#define AB 1
#define ABCW 2
#define ABMM 3
#define ABNN 4

/**************************************************************/
/*                                                            */
/*Procedure de calcul du noyau par appel aux autres procedures*/
/*                                                            */
/**************************************************************/

double kernel(int nbDimensions,individu_t *point1,individu_t *point2,int choixNoyau);

/******************************************************************/
/*                                                                */
/*Procedure de calcul du noyau gaussien                           */
/*bandwidth contient les variances pour les differentes dimensions*/
/*                                                                */
/******************************************************************/

double radialGaussianKernel(int nbDimensions,individu_t *point1,individu_t *point2);


/******************************************************************/
/*                                                                */
/*Procedure de calcul du noyau d'Epanenchnikov                    */
/*bandwidth contient les variances pour les differentes dimensions*/
/*                                                                */
/******************************************************************/

double EpanechnikovKernel(int nbDimensions,individu_t *point1,individu_t *point2);

/******************************************************************/
/*                                                                */
/*Procedure de calcul du noyau polynomial                         */
/*bandwidth contient les variances pour les differentes dimensions*/
/*                                                                */
/******************************************************************/

double polynomialKernel(int nbDimensions,individu_t *point1,individu_t *point2,int choixNoyau);

/***************************************************************/
/*                                                             */
/*Procedure de calcul du bandwidth par estimation de la densite*/
/*                                                             */
/***************************************************************/

void calculFixedBandwidth(int nbDimensions,int nbIndividus,individu_t **individus,double bandwidth);

/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculAdaptiveBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet,
			     int kPlusProches,double scale);

/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculAdaptiveBandwidthComponentWise(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet,
					  int kPlusProches,double scale);

/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculMMDensityBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet) ;

/*bandwidth adaptatif d'Abramson, densite pilote a noyau (kkmeans)*/
void calculPilotAdaptiveBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet);


/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculNNAdaptiveBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet,
			       int kPlusProches,double scale,int accelerationCalculs,int kPlusProchesVoisins); 
