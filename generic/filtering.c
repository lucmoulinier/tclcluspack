#include "mainClusPack.h"
#include "tcl.h"

/******************************************************/
/*                                                    */
/*Procedure de filtrage des dimensions ou celles ayant*/ 
/*des variances sont trop petites sont eliminees      */
/*                                                    */
/******************************************************/

void filtreDimensions(int nbIndividus,
		      int nbDimensionsOrigine,
		      int *nbDimensionsUtilisees,
		      individu_t **individus) 
{
  int i,j;
  double temp,moyenne,variance,varianceLimite=0.00000001;
  
  *nbDimensionsUtilisees=0;
  for(i=0;i<nbDimensionsOrigine;i++) {
    moyenne=0.0;
    for(j=0;j<nbIndividus;j++) {
      moyenne+=individus[j]->valeursBrutes[i];
    }
    moyenne/=(double)nbIndividus;

    variance=0;
    for(j=0;j<nbIndividus;j++) {
      temp=individus[j]->valeursBrutes[i]-moyenne;
      variance+=temp*temp;
    }
    variance/=(double)(nbIndividus-1);
    
    if(variance>varianceLimite) {
      for(j=0;j<nbIndividus;j++) {
        individus[j]->valeursTraitees[*nbDimensionsUtilisees]=individus[j]->valeursBrutes[i];
      }
      (*nbDimensionsUtilisees)++;
    }
  }
}

/******************************************************/
/*                                                    */
/*Procedure de filtrage des dimensions ou celles ayant*/ 
/*des variances sont trop petites sont eliminees      */
/*                                                    */
/******************************************************/

void filtreDimensionsDeuxEnsembles(int nbIndividus,int nbDimensionsOrigine,int *nbDimensionsUtilisees,
				   individu_t **individus,int nbIndividus2,individu_t **individus2) {
  int i,j;
  double temp,moyenne,variance,varianceLimite=0.00000001;

  *nbDimensionsUtilisees=0;
  for(i=0;i<nbDimensionsOrigine;i++) {
    moyenne=0;
    for(j=0;j<nbIndividus;j++) {
      moyenne+=individus[j]->valeursBrutes[i];
    }
    moyenne/=(double)nbIndividus;

    variance=0;
    for(j=0;j<nbIndividus;j++) {
      temp=individus[j]->valeursBrutes[i]-moyenne;
      variance+=temp*temp;
    }
    variance/=(double)(nbIndividus-1);
    
    if(variance>varianceLimite) {
      for(j=0;j<nbIndividus;j++) {
        individus[j]->valeursTraitees[*nbDimensionsUtilisees]=individus[j]->valeursBrutes[i]/sqrt(variance);
      }
      for(j=0;j<nbIndividus2;j++) {
        individus2[j]->valeursTraitees[*nbDimensionsUtilisees]=individus2[j]->valeursBrutes[i]/sqrt(variance);
      }
      (*nbDimensionsUtilisees)++;
    }
  }
}

/*********************************************/
/*                                           */
/*Procedure de filtrage des individus suivant*/ 
/*une valeur seuil et le type des donnees    */
/*                                           */
/*********************************************/

void filtreIndividus(int *nbIndividus,
		     int nbDimensions,
		     individu_t **individus,
		     int typeDonnees,
		     double distanceMin,
		     int distancesDirectes)
{
  /*declaration de variables*/
  int i,j,k,*individusElimines,nbIndividusElimines,compteur;
  double distance;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  individusElimines=(int *)ckalloc(sizeof(int)*(*nbIndividus));
  /*fin allocation memoire*/
  
  nbIndividusElimines=0;
  
  for(i=0;i<*nbIndividus;i++)
    {
      individusElimines[i]=NON;
      individus[i]->nbIndividusSimilaires=0;
      individus[i]->individusSimilaires=NULL;
    }
  for(i=0;i<*nbIndividus;i++)
    {
      for(j=i+1;j<*nbIndividus;j++)
	{
	  if(distancesDirectes==OUI)
	    {
	      /*ALIGNEMENT (bionj) ou DISTANCES : lignes de distances*/
	      distance=individus[i]->valeursTraitees[individus[j]->id];
	    }
	  /*COORDONNEES (y compris alignement projete par MDS)*/
	  else
	    {
	      distance=0;
	      for(k=0;k<nbDimensions;k++)
		{
		  distance+=(individus[i]->valeursTraitees[k]-
			     individus[j]->valeursTraitees[k])*
		    (individus[i]->valeursTraitees[k]-
		     individus[j]->valeursTraitees[k]);
		}
	      distance=sqrt(distance);
	    }
	  if((distance<distanceMin)&&(individusElimines[i]==NON)&&
	     (individusElimines[j]==NON))
	    {
	      nbIndividusElimines++;
	      individusElimines[j]=OUI;
	      individus[i]->nbIndividusSimilaires++;
	      individus[i]->individusSimilaires=(individu_t **)ckrealloc(individus[i]->individusSimilaires,sizeof(individu_t *)*individus[i]->nbIndividusSimilaires);
	      
	      individus[i]->individusSimilaires[individus[i]->nbIndividusSimilaires-1]=individus[j];
	    }
	}
    }

  compteur=0;
  for(i=0;i<*nbIndividus;i++)
    {
      if(individusElimines[i]==NON)
	{
	  individus[compteur]=individus[i];
	  compteur++;
	}
    }
  *nbIndividus-=nbIndividusElimines;

  /*desallocation memoire*/
  ckfree(individusElimines);
  /*fin desallocation memoire*/
}

