#include "mainClusPack.h"
#include "kmeans.h"
#include "tools.h"

/****************************************************/
/*                                                  */
/*Procedure de realisation de clustering par k-means*/
/*en partant d'un debut de clustering               */ 
/*                                                  */
/****************************************************/

int kmeans(int nbIndividus,int nbDimensions,individu_t **individus,int nbClasses)
{
  int i,j,*nbIndividusParClasse,*solution,compteurIterations1,compteurIterations2;
  int nbKmeans,nbIterationsMax,*pointsSelectionnes,resultat;
  double **centres,inertieInterClassePrecedente,inertieInterClasseCourante;
  double meilleureInertieInterClasse,*centreGravite;

  /*allocation memoire*/
  centres=(double **)ckalloc(sizeof(double *)*nbClasses);
  for(i=0;i<nbClasses;i++)
    {
      centres[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  nbIndividusParClasse=(int *)ckalloc(sizeof(int)*nbClasses);
  solution=(int *)ckalloc(sizeof(int)*nbIndividus);
  pointsSelectionnes=(int *)ckalloc(sizeof(int)*nbIndividus);
  centreGravite=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  nbKmeans=3;
  nbIterationsMax=20;
  meilleureInertieInterClasse=NOT_DEFINED;

  /*calcul du centre de gravite de l'ensemble des points*/
  for(i=0;i<nbDimensions;i++) {
    centreGravite[i]=0;
  }
  for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbDimensions;j++) {
      centreGravite[j]+=individus[i]->valeursTraitees[j];
    }
  }
  for(i=0;i<nbDimensions;i++) {
    centreGravite[i]/=(double)nbIndividus;
  }
  
  /*kmeans*/
  for(compteurIterations1=0;compteurIterations1<nbKmeans;compteurIterations1++)
    {
      for(i=0;i<nbClasses;i++)
	{
	  pointsSelectionnes[i]=-1;
	  while(pointsSelectionnes[i]==-1)
	    {
	      pointsSelectionnes[i]=
		(int)floor((double)nbIndividus-1.0)*runif(0,1);
	      for(j=0;j<i;j++)
		{
		  if(pointsSelectionnes[i]==pointsSelectionnes[j])
		    {
		      pointsSelectionnes[i]=-1;
		      break;
		    }
		}
	    }
	  
	  for(j=0;j<nbDimensions;j++)
	    {
	      centres[i][j]=individus[pointsSelectionnes[i]]->valeursTraitees[j];
	    }
	}      
      /*on affecte les individus a leurs groupes les plus proches*/
      resultat=affectationIndividusKmeans(nbDimensions,nbIndividus,nbClasses,nbIndividusParClasse,
					  individus,centres,&inertieInterClassePrecedente,
					  centreGravite);
      
      if(resultat==ECHEC)
	{
	  break;
	}

      compteurIterations2=0;
      while(1)
	{
	  
	  /*on affecte les individus a leurs groupes les plus proches*/
	  resultat=affectationIndividusKmeans(nbDimensions,nbIndividus,nbClasses,
					      nbIndividusParClasse,individus,centres,
					      &inertieInterClasseCourante,centreGravite);
	  
	  if(resultat==ECHEC)
	    {
	      break;
	    }
	  
	  if((-inertieInterClasseCourante+inertieInterClassePrecedente<ZERO_LIMIT)&&
	     (inertieInterClasseCourante-inertieInterClassePrecedente<ZERO_LIMIT))
	    {
	      break;
	    }
	  else
	    {
	      inertieInterClassePrecedente=inertieInterClasseCourante;
	    }

	  compteurIterations2++;
	  
	  if(compteurIterations2>nbIterationsMax)
	    {
	      break;
	    }
	}

      if((resultat!=ECHEC)&&((meilleureInertieInterClasse=NOT_DEFINED)||
	 (inertieInterClasseCourante>meilleureInertieInterClasse)))
	{
	  for(i=0;i<nbIndividus;i++)
	    {
	      solution[i]=individus[i]->cluster;
	    }
	  meilleureInertieInterClasse=inertieInterClasseCourante;
	}
    }
  for(i=0;i<nbIndividus;i++)
    {
      individus[i]->cluster=solution[i];
    }
  
  /*desallocation memoire*/
  for(i=0;i<nbClasses;i++)
    {
      ckfree(centres[i]);
    }
  ckfree(centres);
  ckfree(nbIndividusParClasse);
  ckfree(solution);
  ckfree(pointsSelectionnes);
  ckfree(centreGravite);
  /*fin desallocation memoire*/

  return resultat;
}

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectationIndividusKmeans(int nbDimensions,int nbIndividus,int nbClasses,
				int *nbIndividusParClasse,individu_t **individus,
				double **centres,double *inertieInterClasse,
				double *centreGravite)
{
  /*declaration des variables*/
  int i,j,groupeLePlusProche;
  double distanceCourante,distanceMin,temp1,temp2;
  /*fin declaration des variables*/

  for(i=0;i<nbClasses;i++)
    {
      nbIndividusParClasse[i]=0;
    }

  for(i=0;i<nbIndividus;i++)
    {
      distanceMin=calculDistanceCarre(nbDimensions,centres[0],individus[i]->valeursTraitees);
      groupeLePlusProche=0;

      for(j=1;j<nbClasses;j++)
	{
	  distanceCourante=calculDistanceCarre(nbDimensions,centres[j],
						individus[i]->valeursTraitees);
	  if(distanceCourante<distanceMin)
	    {
	      distanceMin=distanceCourante;
	      groupeLePlusProche=j;
	    }
	}
      individus[i]->cluster=groupeLePlusProche;
      (nbIndividusParClasse[groupeLePlusProche])++;
    }

  for(i=0;i<nbClasses;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  centres[i][j]=0;
	}
    }

  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  centres[individus[i]->cluster][j]+=individus[i]->valeursTraitees[j];
	}
    }
  
  for(i=0;i<nbClasses;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  if(nbIndividusParClasse[i]==0)
	    {
	      return ECHEC;
	    }
	  centres[i][j]/=(double)nbIndividusParClasse[i];
	}
    }
  *inertieInterClasse=0;
  for(i=0;i<nbClasses;i++)
    {
      temp1=0;
      for(j=0;j<nbDimensions;j++)
	{
	  temp2=centres[i][j]-centreGravite[j];
	  temp1+=temp2*temp2;
	}
      *inertieInterClasse+=temp1*(double)nbIndividusParClasse[i];
    }  

  return SUCCES;
}
