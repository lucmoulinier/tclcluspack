#include "mainClusPack.h"
#include "matrice.h"
#include "tcl.h"

/****************************************************/
/*                                                  */
/*Procedure d'estimation de la matrice de covariance*/
/*pour une matrice points x coordonnees             */
/*                                                  */
/****************************************************/

void calculMatriceVariance(int nbDimensions,int nbPoints,double **points,
			   double *matriceVariance)
{
  /*declaration de variables*/
  int i,j,k;
  double *moyennes,diff;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  moyennes=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  /*calcul des moyennes*/
  for(j=0;j<nbDimensions;j++)
    {
      moyennes[j]=0;
    }
  for(i=0;i<nbPoints;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  moyennes[j]+=points[i][j];
	}
    }
  for(j=0;j<nbDimensions;j++)
    {
      moyennes[j]/=(double)nbPoints;
    }

  for(i=0;i<nbDimensions;i++)
    {
      matriceVariance[i]=0;
    }	  
  for(k=0;k<nbPoints;k++)
    {
      for(i=0;i<nbDimensions;i++)
	{
	  diff=points[k][i]-moyennes[i];
	  matriceVariance[i]+=diff*diff;
	}	 
    }
  for(i=0;i<nbDimensions;i++)
    {
      matriceVariance[i]/=(double)(nbPoints-1);
    }	
  
  /*desallocation memoire*/
  ckfree(moyennes);
  /*fin desallocation memoire*/
}

/****************************************************/
/*                                                  */
/*Procedure d'estimation de la matrice de covariance*/
/*                                                  */
/****************************************************/

void calculMatriceCovariance(int nbDimensions,int nbPoints,double **points,
			       double **matriceCovariance)
{
  /*declaration de variables*/
  int i,j,k;
  double *moyennes,diff;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  moyennes=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  /*calcul des moyennes*/
  for(j=0;j<nbDimensions;j++)
    {
      moyennes[j]=0;
    }
  for(i=0;i<nbPoints;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  moyennes[j]+=points[i][j];
	}
    }
  for(j=0;j<nbDimensions;j++)
    {
      moyennes[j]/=(double)nbPoints;
    }

  for(i=0;i<nbDimensions;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  matriceCovariance[i][j]=0;
	}
    }	  
  for(k=0;k<nbPoints;k++)
    {
      for(i=0;i<nbDimensions;i++)
	{
	  diff=points[k][i]-moyennes[i];
	  for(j=i;j<nbDimensions;j++)
	    {
	      matriceCovariance[i][j]+=diff*(points[k][j]-moyennes[j]);
	    }
	}	 
    }
  for(i=0;i<nbDimensions;i++)
    {
      for(j=i;j<nbDimensions;j++)
	{
	  matriceCovariance[i][j]/=(double)(nbPoints-1);
	  
	  matriceCovariance[j][i]=matriceCovariance[i][j];
	}
    }	
  
  
  /*desallocation memoire*/
  ckfree(moyennes);
  /*fin desallocation memoire*/
}


/*******************************************************************/
/*                                                                 */
/*Procedure d'inversion de matrice qui renvoie le determinant      */
/*de la matrice a inverser ou 0 s'il y a un depassement de capacite*/
/*Nota bene : matrice est modifiee par la procedure                */
/*                                                                 */
/*******************************************************************/

double inversionMatrice(int nbDimensions,double **matrice)
{
  int singuliere;
  double logDeterminant;

  logDeterminant=inversionMatriceSinguliere(nbDimensions,matrice,&singuliere);
  if(singuliere) {
    /*comportement historique : 0 si la matrice n'est pas inversible*/
    return 0.0;
  }
  return logDeterminant;
}

/*******************************************************************/
/*                                                                 */
/*Inversion de Gauss-Jordan avec pivot partiel. Renvoie le log du  */
/*determinant ; *singuliere vaut OUI (et matrice est alors         */
/*inutilisable) si un pivot est negligeable devant le plus grand   */
/*element de la matrice.                                           */
/*Nota bene : matrice est modifiee par la procedure                */
/*                                                                 */
/*******************************************************************/

double inversionMatriceSinguliere(int nbDimensions,double **matrice,int *singuliere)
{
  /*declaration de variables*/
  int i,j,k,pivot;
  double **matriceInverse,valeur,echelle;
  double determinant;
  /*fin declaration de variables*/

  *singuliere=NON;

  /*allocation memoire*/
  matriceInverse=(double **)ckalloc(sizeof(double *)*nbDimensions);
  for(i=0;i<nbDimensions;i++)
    {
      matriceInverse[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  /*fin allocation memoire*/

  echelle=0.0;
  for(i=0;i<nbDimensions;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  matriceInverse[i][j]=0;
	  if(fabs(matrice[i][j])>echelle) {
	    echelle=fabs(matrice[i][j]);
	  }
	}
      matriceInverse[i][i]=1.0;
    }

  determinant=0.0;
  for(i=0;i<nbDimensions;i++)
    { 
      /*pivot partiel : plus grand element de la colonne*/
      pivot=i;
      for(j=i+1;j<nbDimensions;j++)
	{
	  if(fabs(matrice[j][i])>fabs(matrice[pivot][i]))
	    {
	      pivot=j;
	    }
	}
      
      if(!(fabs(matrice[pivot][i])>ZERO_LIMIT*echelle))
	{
	  *singuliere=OUI;
	  determinant=-HUGE_VAL;
	  goto fin;
	}
      echangeLignes(nbDimensions,matrice,i,pivot);
      echangeLignes(nbDimensions,matriceInverse,i,pivot);
 
      valeur=matrice[i][i];

      for(j=0;j<nbDimensions;j++)
	{
	  matrice[i][j]/=valeur;
	  matriceInverse[i][j]/=valeur;
	}
      determinant+=log(fabs(valeur));
      
      for(j=0;j<nbDimensions;j++)
	{
	  if(j!=i)
	    { 
	      valeur=matrice[j][i];
	      for(k=0;k<nbDimensions;k++)
		{
		  matrice[j][k]-=valeur*matrice[i][k];
		  matriceInverse[j][k]-=valeur*matriceInverse[i][k];
		}
	    }	
	}

    }
  
  for(i=0;i<nbDimensions;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  matrice[i][j]=matriceInverse[i][j];
	}
    }

 fin:
  /*desallocation memoire*/
  for(i=0;i<nbDimensions;i++)
    {
      ckfree((char *)matriceInverse[i]);
    }
  ckfree((char *)matriceInverse);
  /*fin desallocation memoire*/
  
  return determinant;
}

/*******************************************************************/
/*                                                                 */
/*Inverse une matrice de covariance source (non modifiee) dans     */
/*inverse et renvoie le log du determinant. Si la matrice est      */
/*singuliere (ex. : donnees de composition dont la somme vaut 1),  */
/*on ajoute a la diagonale un terme lambda (1e-6 x variance        */
/*moyenne, multiplie par 10 a chaque echec).                       */
/*                                                                 */
/*******************************************************************/

double inversionMatriceRegularisee(int nbDimensions,double **source,double **inverse)
{
  int i,j,essai,singuliere;
  double logDeterminant,lambda,trace;

  trace=0.0;
  for(i=0;i<nbDimensions;i++) {
    trace+=fabs(source[i][i]);
  }
  lambda=(trace>0.0) ? 1e-6*trace/(double)nbDimensions : 1e-6;

  for(essai=0;;essai++) {
    for(i=0;i<nbDimensions;i++) {
      for(j=0;j<nbDimensions;j++) {
	inverse[i][j]=source[i][j];
      }
      if(essai>0) {
	inverse[i][i]+=lambda;
      }
    }
    logDeterminant=inversionMatriceSinguliere(nbDimensions,inverse,&singuliere);
    if(!singuliere||essai>=12) {
      break;
    }
    if(essai>0) {
      lambda*=10.0;
    }
  }
  if(singuliere) {
    /*ne devrait pas arriver : matrice identite*/
    for(i=0;i<nbDimensions;i++) {
      for(j=0;j<nbDimensions;j++) {
	inverse[i][j]=(i==j) ? 1.0 : 0.0;
      }
    }
    logDeterminant=0.0;
  }
  return logDeterminant;
}

/*****************************************************/
/*                                                   */
/*Procedure d'echange de deux lignes dans une matrice*/
/*                                                   */
/*****************************************************/

void echangeLignes(int nbDimensions,double **matrice,int ligne1,int ligne2)
{
  /*declaration de variables*/
  int i;
  double valTemp;
  /*fin declaration de variables*/

  for(i=0;i<nbDimensions;i++)
    {
      valTemp=matrice[ligne1][i];
      matrice[ligne1][i]=matrice[ligne2][i];
      matrice[ligne2][i]=valTemp;
    }
}


/*******************************************************************************/
/*                                                                             */
/*Procedure de calcul de la norme au carre d'un vecteur pour une matrice donnee*/
/*                                                                             */
/*******************************************************************************/

double calculNorme(int nbDimensions,double *vecteur,double **matrice)
{
  /*declaration de variables*/
  int i,j;
  double normeCarre;
  /*fin declaration de variables*/
  
  /*calcul du produit de gauche*/
  normeCarre=0;
  for(i=0;i<nbDimensions;i++)
    {
      for(j=i+1;j<nbDimensions;j++)
	{
	  normeCarre+=2.0*vecteur[i]*vecteur[j]*matrice[j][i];
	}
      normeCarre+=vecteur[i]*vecteur[i]*matrice[i][i];
    }

  return normeCarre;
}

/******************************************/
/*                                        */
/*Procedure effectuant vecteur2=A*vecteur1*/
/*                                        */
/******************************************/

void multiplicationMatriceVecteur(int nbDimensions,double **A,double *vecteur1,double *vecteur2) {
  int i,j;

  for(i=0;i<nbDimensions;i++) {
    vecteur2[i]=0;
    for(j=0;j<nbDimensions;j++) {
      vecteur2[i]+=A[i][j]*vecteur1[j];
    }
  }
}

/**************************************************************************/
/*                                                                        */
/*Procedure calculant le vecteur propre associe a une valeur propre donnee*/
/*                                                                        */
/**************************************************************************/

void calculVecteurPropre(int nbDimensions,double **matrice,double lambda,double *vecteurPropre) {
  int i,j;
  double **A,**B;

  /*allocation memoire*/
  A=(double **)ckalloc(sizeof(double *)*nbDimensions);
  B=(double **)ckalloc(sizeof(double *)*nbDimensions);
  for(i=0;i<nbDimensions;i++) {
    A[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    B[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }  
  /*fin allocation memoire*/

  for(i=0;i<nbDimensions;i++) {
    for(j=0;j<nbDimensions;j++) {
      A[i][j]=matrice[i][j];
    }
    A[i][i]-=lambda;
  }

  /*on cherche une variable pouvant etre fixee a 1,
   pour cela on cherche une variable qui n'est jamais seule dans une equation*/
  int variableJamaisSeule,queDesZeros,variable;
  for(variable=0;variable<nbDimensions;variable++) {
    variableJamaisSeule=OUI;
    for(i=0;i<nbDimensions;i++) { 
      queDesZeros=OUI;
      for(j=0;j<nbDimensions;j++)
	{
	  if(j!=variable) {
	    if(fabs(A[i][j])>ZERO_LIMIT)
	      {
		queDesZeros=NON;
		break;
	      }
	  }
	}
      if(queDesZeros==OUI) {
	variableJamaisSeule=NON;
	break;
      }
    }  
    if(variableJamaisSeule==OUI) {
      break;
    }
  }
  for(j=0;j<variable;j++) {
    for(i=0;i<nbDimensions-1;i++) {
      B[i][j]=A[i][j];
    }
  }
  for(j=variable+1;j<nbDimensions;j++) {
    for(i=0;i<nbDimensions-1;i++) {
      B[i][j-1]=A[i][j];
    }
  }

  /*inversion de la matrice*/
  inversionMatrice(nbDimensions-1,B);

  /*calcul du vecteur propre*/
  for(i=0;i<variable;i++) {
    vecteurPropre[i]=0;
    for(j=0;j<nbDimensions-1;j++) {
      vecteurPropre[i]+=-B[i][j]*A[j][variable];
    }
  }
  vecteurPropre[variable]=1.0;
  for(i=variable+1;i<nbDimensions;i++) {
    vecteurPropre[i]=0;
    for(j=0;j<nbDimensions-1;j++) {
      vecteurPropre[i]+=-B[i-1][j]*A[j][variable];
    }
  }
  
  /*desallocation memoire*/
  for(i=0;i<nbDimensions;i++) {
    ckfree(A[i]);
    ckfree(B[i]);
  }
  ckfree(A);
  ckfree(B);
  /*fin desallocation memoire*/
}
