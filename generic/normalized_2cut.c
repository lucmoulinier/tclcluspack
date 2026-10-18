#include "mainClusPack.h"
#include "symeig.h"
#include "tools.h"
#include "normalized_2cut.h"

/********************************************************************/
/*                                                                  */
/*Procedure de recherche de coupe minimum normalisee de Xing et Karp*/
/*                                                                  */
/********************************************************************/

int normalizedCut(cluster_t *groupe,cluster_t *groupe1,cluster_t *groupe2)
{
  /*declaration de variables*/
  int i,j,nbIndividus,*positions,meilleureSolution;
  double **DminusW,**A,*D,*W,**V,*vecteurPropre,meilleurScore,score,*y;
  double *coordonneesTriees;
  /*fin declaration de variables*/

  nbIndividus=groupe->nbIndividus;

  /*allocation memoire*/
  DminusW=(double **)ckalloc(sizeof(double *)*nbIndividus);
  A=(double **)ckalloc(sizeof(double *)*nbIndividus);
  V=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      DminusW[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
      A[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
      V[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
    }
  D=(double *)ckalloc(sizeof(double)*nbIndividus);
  W=(double *)ckalloc(sizeof(double)*nbIndividus);
  y=(double *)ckalloc(sizeof(double)*nbIndividus);
  vecteurPropre=(double *)ckalloc(sizeof(double)*nbIndividus);
  coordonneesTriees=(double *)ckalloc(sizeof(double)*nbIndividus);
  positions=(int *)ckalloc(sizeof(int)*nbIndividus);
  /*fin allocation memoire*/

  for(i=0;i<nbIndividus;i++)
    {
      D[i]=0;
      for(j=0;j<nbIndividus;j++)
	{
	  D[i]+=(double)groupe->individus[i]->valeursTraitees[groupe->individus[j]->id];
	  DminusW[i][j]=0;
	}
    }

  for(i=0;i<nbIndividus;i++)
    {
      DminusW[i][i]=D[i];
      for(j=i;j<nbIndividus;j++)
	{
	  DminusW[i][j]-=(double)groupe->individus[i]->valeursTraitees[groupe->individus[j]->id];
	  A[i][j]=DminusW[i][j];
	  DminusW[j][i]=DminusW[i][j];
	  A[j][i]=A[i][j];
	}
    }

  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbIndividus;j++)
	{
	  A[i][j]/=sqrt(D[i]);
	}
    }

  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbIndividus;j++)
	{
	  A[j][i]/=sqrt(D[i]);
	}
    }

  /*on cherche les vecteurs propres*/
  symeig(A,nbIndividus,V,W);

  /*on cherche le vecteur propre a la deuxieme plus petite valeur propre*/
  for(i=0;i<nbIndividus;i++)
    {
      positions[i]=i;
    }

  /*on trie dans l'ordre decroissant les valeurs propres*/
  triRapideValeursPositions(W,positions,0,nbIndividus-1);
  
  for(i=0;i<nbIndividus;i++)
    {
      vecteurPropre[i]=V[i][positions[nbIndividus-2]]/sqrt(D[i]);
      coordonneesTriees[i]=vecteurPropre[i];
    }

  for(i=0;i<nbIndividus;i++)
    {
      positions[i]=i;
    }

  /*on trie dans l'ordre decroissant les coordonnees du vecteur propre*/
  triRapideValeursPositions(coordonneesTriees,positions,0,nbIndividus-1);

  for(i=0;i<nbIndividus;i++)
    {
      y[i]=-1.0;
    }

  meilleurScore=-VALEUR_ENORME;
  for(i=0;i<nbIndividus;i++)
    {
      y[positions[i]]=1.0;
      score=evaluationCoupeNormaliseeMinimum(nbIndividus,y,D,DminusW);
      printf("score : %f\n",score);
      if(score>meilleurScore)
	{
	  meilleurScore=score;
	  meilleureSolution=i;
	  printf("i : %d \n",i);
	}
    }

  for(i=0;i<nbIndividus;i++)
    {
      y[i]=-1.0;
    }
  for(i=0;i<=meilleureSolution;i++)
    {
      y[positions[i]]=1.0;
    }
  
  groupe1->nbIndividus=0;
  groupe2->nbIndividus=0;
  for(i=0;i<nbIndividus;i++)
    {
      if(y[i]<0)
	{
	  groupe1->individus[groupe1->nbIndividus]=groupe->individus[i];
	  groupe1->nbIndividus++;
	}
      else
	{
	  groupe2->individus[groupe2->nbIndividus]=groupe->individus[i];
	  groupe2->nbIndividus++;
	}
    }
  printf("groupe1 : %d ; groupe2 : %d\n",groupe1->nbIndividus,groupe2->nbIndividus);

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++)
    {
      ckfree(DminusW[i]);
      ckfree(A[i]);
      ckfree(V[i]);
    }
  ckfree(DminusW);
  ckfree(A);
  ckfree(D);
  ckfree(V);
  ckfree(W);
  ckfree(vecteurPropre);
  ckfree(y);
  ckfree(coordonneesTriees);
  ckfree(positions);
  /*fin desallocation memoire*/

  return OUI;
}

/********************************************************/
/*                                                      */
/*Procedure renvoyant pour un vecteur donne la valeur   */
/*de la fonction objectif associee a la coupe normalisee*/
/*                                                      */
/********************************************************/

double evaluationCoupeNormaliseeMinimum(int nbDimensions,double *y,double *D,double **DminusW)
{
  /*declaration de variables*/
  int i,j;
  double res,numerateur,denominateur,*vecteur;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  vecteur=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  denominateur=0;
  for(i=0;i<nbDimensions;i++)
    {
      denominateur+=y[i]*y[i]*D[i];
    }

  for(i=0;i<nbDimensions;i++)
    {
      vecteur[i]=0;
      for(j=0;j<nbDimensions;j++)
	{
	  vecteur[i]+=y[j]*DminusW[j][i];
	}
    }

  numerateur=0;
  for(i=0;i<nbDimensions;i++)
    {
      numerateur+=vecteur[i]*y[i];
    }

  res=numerateur/denominateur;

  /*desallocation memoire*/
  ckfree(vecteur);
  /*fin desallocation memoire*/

  return res;
}

