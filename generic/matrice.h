/*******************************************************************************/
/*                                                                             */
/*Procedure de calcul de la norme au carre d'un vecteur pour une matrice donnee*/
/*                                                                             */
/*******************************************************************************/

double calculNorme(int nbDimensions,double *vecteur,double **matrice);

/****************************************************/
/*                                                  */
/*Procedure d'estimation de la matrice de covariance*/
/*pour une matrice points x coordonnees             */
/*                                                  */
/****************************************************/

void calculMatriceVariance(int nbDimensions,int nbPoints,double **points,
			     double *matriceVariance);

/****************************************************/
/*                                                  */
/*Procedure d'estimation de la matrice de covariance*/
/*                                                  */
/****************************************************/

void calculMatriceCovariance(int nbDimensions,int nbPoints,double **points,
			       double **matriceCovariance);

/*******************************************************************/
/*                                                                 */
/*Procedure d'inversion de matrice qui renvoie le determinant      */
/*de la matrice a inverser ou 0 s'il y a un depassement de capacite*/
/*Nota bene : matrice est modifiee par la procedure                */
/*                                                                 */
/*******************************************************************/

double inversionMatrice(int nbDimensions,double **matrice);
double inversionMatriceSinguliere(int nbDimensions,double **matrice,int *singuliere);
double inversionMatriceRegularisee(int nbDimensions,double **source,double **inverse);

/*****************************************************/
/*                                                   */
/*Procedure d'echange de deux lignes dans une matrice*/
/*                                                   */
/*****************************************************/

void echangeLignes(int nbDimensions,double **matrice,int ligne1,int ligne2);


/*****************************************/
/*                                       */
/*Procedure effectuan vecteur2=A*vecteur1*/
/*                                       */
/*****************************************/

void multiplicationMatriceVecteur(int nbDimensions,double **A,double *vecteur1,double *vecteur2);

/**************************************************************************/
/*                                                                        */
/*Procedure calculant le vecteur propre associe a une valeur propre donnee*/
/*                                                                        */
/**************************************************************************/

void calculVecteurPropre(int nbDimensions,double **matrice,double lambda,double *vecteurPropre);
