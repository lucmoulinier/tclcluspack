#include "mainClusPack.h"
#include "tools.h"

/**************************************************************/
/*                                                            */
/*Procedure de calcul du volume d'un tetraedre a nb_dimensions*/
/*                                                            */
/**************************************************************/

void calculVolumeSimplex(int nbDimensions,double *volumeSimplex)
{
  /*declaration de variables*/
  int i;
  /*fin declaration de variables*/

  *volumeSimplex=1.0;
  for(i=1;i<=nbDimensions;i++)
    {
      *volumeSimplex/=(double)i;
    }
  *volumeSimplex*=sqrt((double)(nbDimensions+1));
}

/*************************************************************************************/
/*                                                                                   */
/*Procedure de calcul d'une sphere unitaire a nbDimensions, dans R^{nbDimensions+1}*/
/*                                                                                   */
/*************************************************************************************/

void calculSurfaceSphere(int nbDimensions,double *surfaceSphere)
{
  /*declaration de variables*/
  int i;
  /*fin declaration de variables*/

  if(nbDimensions%2==0)
    {
      *surfaceSphere=2.0;
      for(i=1;i<=nbDimensions/2;i++)
	{
	  (*surfaceSphere)*=2.0*M_PI/(double)(2*i-1);
	}
    }
  else
    {
      *surfaceSphere=2.0*M_PI;
      for(i=1;i<=(nbDimensions-1)/2;i++)
	{
	  (*surfaceSphere)*=M_PI/(double)i;
	}
    }
}


/*********************************************************************************************/
/*                                                                                           */
/*Procedure de calcul de comparaison pour le classement des individus par distance croissante*/
/*                                                                                           */
/*********************************************************************************************/

int compareIndividus(void const *a,void const *b) {
  individu_t **pa=(individu_t **)a;
  individu_t **pb=(individu_t **)b;
  double diff=(*pa)->valeurTemp-(*pb)->valeurTemp;

  if(diff>0) {
    return 1;
  } else {
    return 0;
  }
}

/********************************************************************************************/
/*                                                                                          */
/*Procedure de calcul de la densite non parametrique en utilisant les k plus proches voisins*/
/*                                                                                          */
/********************************************************************************************/

double calculDensiteNonParametrique(cluster_t *groupe,double *coordonnees,int dimension,
				    int nbVoisins,int ajout) {
  /*declaration de variables*/
  int i,j,borneSup;
  individu_t *individuCourant;
  double densite,volume,min,max;
  /*fin declaration de variables*/

  for(i=0;i<groupe->nbIndividus;i++) {
    individuCourant=groupe->individus[i];
    individuCourant->valeurTemp=
      calculDistanceCarre(dimension,coordonnees,individuCourant->valeursTraitees);
  }

  /*on trie par ordre croissant*/
  qsort(groupe->individus,groupe->nbIndividus,sizeof(individu_t *),compareIndividus);
  /*  for(i=0;i<groupe->nbIndividus;i++) {
    printf("val : %f\n",groupe->individus[i]->valeurTemp);
  getchar();
  }*/
  volume=1.0;
  for(i=0;i<dimension;i++) {
    min=coordonnees[i];
    max=coordonnees[i];
    borneSup=nbVoisins;
    if(borneSup>groupe->nbIndividus) {
      borneSup=groupe->nbIndividus;
    }
    for(j=0;j<borneSup;j++) {
      if(groupe->individus[j]->valeursTraitees[i]<min) {
	min=groupe->individus[j]->valeursTraitees[i];
      } else if(groupe->individus[j]->valeursTraitees[i]>max) {
	max=groupe->individus[j]->valeursTraitees[i];
      }
    }
    volume*=(max-min);
  }
  densite=(double)(nbVoisins+ajout)/volume;


  return densite;
}


/*************************************************************/
/*                                                           */
/*Procedure de calcul des exposants d'Holder                 */
/*Methode : on forme un hyperpave avec pour longueur         */
/*la distances entre les deux coordonnees                    */
/*les plus eloignees sur chaque axe                          */
/*                                                           */
/*************************************************************/

void calculExposantsHolderPoints(int nbDimensions,int nbIndividusTotal,
				 double **coordonneesTousIndividus,int nbIndividus,
				 double **coordonneesIndividus,double *exposantsHolder,
				 double *cotesHyperpave,double *valeursMinimum,
				 int modeVerbose)
{
  /*declaration de variables*/
  int i,j,k,l,testIndividuDansLeMemeHyperpave,nbMesures,nbMesuresMax;
  int *individusATester,nbIndividusATester,compteur;
  double *bornesInf,*bornesSup,indiceSegmentAxe,nbCopains;
  /*fin declaration de variables*/

  nbMesuresMax=20;

  /*allocation memoire*/
  bornesInf=(double *)ckalloc(sizeof(double)*nbDimensions);
  bornesSup=(double *)ckalloc(sizeof(double)*nbDimensions);
  individusATester=(int *)ckalloc(sizeof(int)*nbIndividusTotal);
  /*fin allocation memoire*/
  
  for(i=0;i<nbIndividus;i++)
    {
      nbMesures=1;
      for(j=0;j<nbIndividusTotal;j++)
	{
	  individusATester[j]=j;
	}
      nbIndividusATester=nbIndividusTotal;
      for(j=0;j<nbMesuresMax;j++)      
	{
	  for(k=0;k<nbDimensions;k++)
	    {
	      indiceSegmentAxe=floor((double)(coordonneesIndividus[i][k]-valeursMinimum[k])/
				       ((double)cotesHyperpave[k]/pow(2.0,(double)(j+1))));
	   if(indiceSegmentAxe>pow(2.0,(double)(j+1))-ZERO_LIMIT)
		{
		  indiceSegmentAxe-=1.0;
		}
	      bornesInf[k]=indiceSegmentAxe*cotesHyperpave[k]/pow(2.0,(double)(j+1))+
		valeursMinimum[k];
	      bornesSup[k]=(1.0+indiceSegmentAxe)*cotesHyperpave[k]/pow(2.0,(double)(j+1))+
		valeursMinimum[k];   
	    }
	  nbCopains=0;

	  compteur=0;

	  for(k=0;k<nbIndividusATester;k++)
	    {
	      testIndividuDansLeMemeHyperpave=OUI;
	      for(l=0;l<nbDimensions;l++)
		{
		  if((coordonneesTousIndividus[individusATester[k]][l]<
		      bornesInf[l]-ZERO_LIMIT)||
		     (coordonneesTousIndividus[individusATester[k]][l]>
		      bornesSup[l]+ZERO_LIMIT))
		    {
		      testIndividuDansLeMemeHyperpave=NON;
		      break;
		    }
		}
	      if(testIndividuDansLeMemeHyperpave==OUI)
		{
		  nbCopains+=1.0;
		  if(nbCopains>=2.0)
		    {
		      for(l=k;l<nbIndividusATester;l++)
			{
			  individusATester[compteur]=individusATester[l];
			  compteur++;
			}
		      break;
		    }
		  else
		    {
		      individusATester[compteur]=individusATester[k];
		      compteur++;
		    }
		}
	    }
	  nbIndividusATester=compteur;

	  if(modeVerbose==OUI)
	    {
	      printf("nb_copains : %.2f\n",nbCopains);
	    }	 

	  if(nbCopains<1.0)
	    {
	      break;
	    }
	  else if(nbCopains<2.0)
	    {
	      nbMesures++;
	      break;
	    }
	  else
	    {
	      nbMesures++;
	    }
	}

      /*      printf("nbMesures : %d\n",nbMesures);*/
      exposantsHolder[i]=(double)nbMesures;
     
    }

  /*desallocation memoire*/
  ckfree(bornesInf);
  ckfree(bornesSup);
  ckfree(individusATester);
  /*fin desallocation memoire*/
}


