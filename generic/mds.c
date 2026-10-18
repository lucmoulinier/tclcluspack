#include "mainClusPack.h"
#include "mds.h"
#include "svd0.h"
#include "symeig.h"

/***************************************************************************************/
/*                                                                                     */
/*Procedure de calcul des matrices A et B obtenues a partir de la matrice des individus*/
/*                                                                                     */
/***************************************************************************************/     

void matriceAB(int nbIndividus,double **individus,double **A,double **B)
{
  /*declaration des variables*/
  int i,j;
  double *Ai_bar,Aij_bar;
  /*fin de declaration des variables*/

  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbIndividus;j++)
	{
	  A[i][j]=-(individus[i][j]*individus[i][j])/2.0;
	} 
    } 

  /*allocation de memoire pour Ai_bar*/
  Ai_bar=(double *)ckalloc(sizeof(double)*nbIndividus);
  /*fin d'allocation de memoire*/
  
  /*calcul de la moyenne des elements de A sur chaque colonne*/  
  for(i=0;i<nbIndividus;i++)
    {
      Ai_bar[i]=0;
      for(j=0;j<nbIndividus;j++)
	{
	  Ai_bar[i]+=A[j][i];
	}
      Ai_bar[i]/=(double)nbIndividus;
    }

  Aij_bar=0;
  /*calcul de la moyenne de tous les elts de A*/
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbIndividus;j++)
	{
	  Aij_bar+=A[i][j];
	}
    }
  Aij_bar/=(double)(nbIndividus*nbIndividus);
  
  /*calcul des elts de B*/
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbIndividus;j++)
	{
	  B[i][j]=A[i][j]-Ai_bar[i]-Ai_bar[j]+Aij_bar;
	}
    }

  /*desallocation memoire*/
  ckfree(Ai_bar);
  /*fin desallocation memoire*/
}   

/*********************************************************************/
/*                                                                   */
/*Procedure de calcul des coordonnees normalisees qui viennent de MDS*/
/*                                                                   */
/*********************************************************************/

void distancesToCoordonnees(int nbIndividus,
			    individu_t **individus,
			    double **distances,
			    int nbDimensions) 
{
  /*declaration des variables*/
  int i,j;
  double **U,**V,**A,**B,*W;
  /*fin declaration des variables*/
  
  /*allocation de  memoire pour U,V et W de svd0 pour la matrice B de MDS*/
  U=(double **)ckalloc(sizeof(double*)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      U[i]=(double*)ckalloc(sizeof(double)*nbIndividus);
    }
  V=(double **)ckalloc(sizeof(double*)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      V[i]=(double*)ckalloc(sizeof(double)*nbIndividus);
    }
  W=(double *)ckalloc(sizeof(double)*nbIndividus);
  A=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    { 
      A[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
    }
  B=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    { 
      B[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
    }
  /*fin d'allocation de memoire*/
  
  /*calcul des matrices A et B a partir de la matrice de individus*/
  matriceAB(nbIndividus,distances,A,B);

  /*appel de la procedure svd*/
  symeig(B,nbIndividus,U,W);

  /*calcul des coordonnees MDS*/ 
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  individus[i]->valeursTraitees[j]=U[i][nbIndividus-1-j];
	  individus[i]->valeursBrutes[j]=individus[i]->valeursTraitees[j];
	} 
    }

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++)
    {
      ckfree(U[i]);
      ckfree(V[i]);
      ckfree(A[i]);
      ckfree(B[i]);
    }
  ckfree(U);
  ckfree(V);
  ckfree(A);
  ckfree(B);
  ckfree(W);
  /*fin desallocation memoire*/
}
