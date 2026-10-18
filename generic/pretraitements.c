#include "mainClusPack.h"
#include "emSphericalGaussian.h"
#include "symeig.h"
#include "kernels.h"
#include "matrice.h"
#include "dissimilarities_clustering.h"
#include "tcl.h"

/********************************************/
/*                                          */
/*Procedure de standardisation des individus*/
/*                                          */
/********************************************/

int standardisationIndividus(int nbIndividus,int nbDimensions,individu_t **individus)
{
  /*declaration de variables*/
  int i,j;
  double moyenne,ecartType;
  /*fin declaration de variables*/
  
  for(i=0;i<nbIndividus;i++)
    {
      moyenne=0;
      for(j=0;j<nbDimensions;j++)
	{
	  moyenne+=individus[i]->valeursTraitees[j];
	}
      moyenne/=(double)nbDimensions;

      ecartType=0;
      for(j=0;j<nbDimensions;j++)
	{
	  ecartType+=(individus[i]->valeursTraitees[j]-moyenne)*(individus[i]->valeursTraitees[j]-moyenne);
	}
      ecartType=sqrt(ecartType/(double)nbDimensions);

      if(ecartType<ZERO_LIMIT)
	{
	  /*pas d'exit() dans une librairie : l'appelant signale l'erreur*/
	  return i;
	}

      for(j=0;j<nbDimensions;j++)
	{
	  individus[i]->valeursTraitees[j]=(individus[i]->valeursTraitees[j]-moyenne)/ecartType;
	  individus[i]->valeursBrutes[j]=individus[i]->valeursTraitees[j];
	}
    }
  return -1;
}

/*************************************/
/*                                   */
/*Procedure de reduction de dimension*/
/*                                   */
/*************************************/

void reductionDimension(int nbIndividus,individu_t **individus,int nbDimensions,
			int nbDimensionsDemandees) {
  /*declaration de variables*/
  int i,j,k;
  double temp,*moyenne,**matriceCovariance,**vecteursPropres,*valeursPropres;;
  double **nouvellesCoordonnees;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  moyenne=(double *)ckalloc(sizeof(double)*nbDimensions);
  matriceCovariance=(double **)ckalloc(sizeof(double *)*nbDimensions);
  for(i=0;i<nbDimensions;i++) {
    matriceCovariance[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  valeursPropres=(double *)ckalloc(sizeof(double)*nbDimensions);
  vecteursPropres=(double **)ckalloc(sizeof(double *)*nbDimensions);
  for(i=0;i<nbDimensions;i++) {
    vecteursPropres[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  nouvellesCoordonnees=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    nouvellesCoordonnees[i]=(double *)ckalloc(sizeof(double)*nbDimensionsDemandees);
  }
  /*fin allocation memoire*/
  
  /*calcul de la moyenne*/
  for(j=0;j<nbDimensions;j++)
    {
      moyenne[j]=0;
    }
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  moyenne[j]+=individus[i]->valeursTraitees[j];
	}
    }
  for(j=0;j<nbDimensions;j++)
    {
      moyenne[j]/=(double)nbIndividus;
    }
  
  /*calcul de la matrice de variance-covariance*/ 
  for(i=0;i<nbDimensions;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  matriceCovariance[i][j]=0;
	}
    }
  for(k=0;k<nbIndividus;k++)
    {
      for(i=0;i<nbDimensions;i++)
	{ 
	  temp=individus[k]->valeursTraitees[i]-moyenne[i];
	  for(j=i;j<nbDimensions;j++)
	    {
	      matriceCovariance[i][j]+=temp*(individus[k]->valeursTraitees[j]-moyenne[j]);
	    }
	}	 
    }
   for(i=0;i<nbDimensions;i++)
    {
      for(j=i;j<nbDimensions;j++)
	{
	  matriceCovariance[i][j]/=(double)(nbIndividus-1);
	  
	  matriceCovariance[j][i]=matriceCovariance[i][j];
	}
    }	

   /*Principal Component Analysis*/
   /*on calcule les vecteurs propres de la matrice de covariance*/
   symeig(matriceCovariance,nbDimensions,vecteursPropres,valeursPropres);

   for(i=0;i<nbDimensions;i++) {
     printf("eigenvalues : %f\n",valeursPropres[i]);
   }

   /*calcul des nouvelles coordonnees*/
   for(i=0;i<nbIndividus;i++) {
     for(j=0;j<nbDimensionsDemandees;j++) {
       nouvellesCoordonnees[i][j]=0;
       for(k=0;k<nbDimensions;k++) {
	 nouvellesCoordonnees[i][j]+=vecteursPropres[k][nbDimensions-1-j]*
	   (individus[i]->valeursTraitees[k]-moyenne[k]);
       }
       nouvellesCoordonnees[i][j]+=moyenne[j];
     }
   }

   /*mise a jour des coordonnees des individus*/
   for(i=0;i<nbIndividus;i++) {
     for(j=0;j<nbDimensionsDemandees;j++) {
       individus[i]->valeursTraitees[j]=nouvellesCoordonnees[i][j];
     }
   }

   /*desallocation memoire*/
   ckfree(moyenne);
   ckfree(valeursPropres);
   for(i=0;i<nbDimensions;i++) {
     ckfree(matriceCovariance[i]);
     ckfree(vecteursPropres[i]);
   }
   ckfree(matriceCovariance);
   ckfree(vecteursPropres);
   for(i=0;i<nbIndividus;i++) {
     ckfree(nouvellesCoordonnees[i]);
   }
   ckfree(nouvellesCoordonnees);
   /*fin desallocation memoire*/
}

/**************************************************/
/*                                                */
/*Procedure de projection des points sur la sphere*/
/*                                                */
/**************************************************/

void projectionSurSphere(int nbIndividus,int nbDimensions,individu_t **individus) {
  /*declaration de variables*/
  int i,j;
  double rayon,*angles,*vecteurRotation,**nouvellesCoordonnees,**coordonneesIndividus;
  /*fin declaration de variables*/

  /*allocation memoire*/
  angles=(double *)ckalloc(sizeof(double)*(nbDimensions-1));
  vecteurRotation=(double *)ckalloc(sizeof(double)*nbDimensions);
  nouvellesCoordonnees=(double **)ckalloc(sizeof(double *)*nbIndividus);
  coordonneesIndividus=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      nouvellesCoordonnees[i]=(double *)ckalloc(sizeof(double)*nbDimensions); 
      coordonneesIndividus[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  /*fin allocation memoire*/

  /*copie des coordonnees d'individus vers coordonnees individus*/
  for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbDimensions;j++) {
      coordonneesIndividus[i][j]=individus[i]->valeursTraitees[j];
    }
  }

  /*calcul des nouvelles coordonnees*/
  for(i=0;i<nbDimensions;i++)
    {
      vecteurRotation[i]=1.0/sqrt((double)nbDimensions);
    }
  cartesiennes2polaires(nbDimensions,vecteurRotation,angles);
  
  calculCoordonneesApresRotation(nbDimensions,nbIndividus,coordonneesIndividus,
				    angles,nouvellesCoordonnees);
  
  /*on ramene les points sur une sphere de rayon 1*/
  for(i=0;i<nbIndividus;i++)
    {
      rayon=0;
      for(j=0;j<nbDimensions;j++)
	{
	  rayon+=pow(nouvellesCoordonnees[i][j],2.0);
	}
      rayon=sqrt(rayon);
     for(j=0;j<nbDimensions;j++)
	{
	  nouvellesCoordonnees[i][j]/=rayon;
	}
    }

  /*copie de nouvellesCoordonnees dans les coordonnees d'individus*/
  for(i=0;i<nbIndividus;i++) {
    //printf("%s : ",individus[i]->nom);
    for(j=0;j<nbDimensions-1;j++) {
      individus[i]->valeursTraitees[j]=nouvellesCoordonnees[i][j+1];
    }
    //    printf("\n");
  }

  /*desallocation memoire*/
  ckfree(angles);
  ckfree(vecteurRotation);
  for(i=0;i<nbIndividus;i++) {
    ckfree(nouvellesCoordonnees[i]);
    ckfree(coordonneesIndividus[i]);
  }
  ckfree(nouvellesCoordonnees);
  ckfree(coordonneesIndividus);
  /*fin desallocation memoire*/
}

/**************************************************************/
/*                                                            */
/*Procedure de calcul des similarites a partir des coordonnees*/
/*en utilisant les correlations                               */
/*                                                            */
/**************************************************************/

void coordonneesToSimilarites(int nbIndividus,int nbDimensions,individu_t **individus)
{
  /*declaration de variables*/
  int i,j,k;
  double **correlations,*ecartsType,*moyennes;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  correlations=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      correlations[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
    }
  moyennes=(double *)ckalloc(sizeof(double)*nbIndividus);
  ecartsType=(double *)ckalloc(sizeof(double)*nbIndividus);
  /*fin allocation memoire*/
  
  for(i=0;i<nbIndividus;i++)
    {
      moyennes[i]=0;
      for(j=0;j<nbDimensions;j++)
	{
	  moyennes[i]+=individus[i]->valeursTraitees[j];
	}
      moyennes[i]/=(double)nbDimensions;
    }

  for(i=0;i<nbIndividus;i++)
    {
      ecartsType[i]=0;
      for(j=0;j<nbDimensions;j++)
	{
	  ecartsType[i]+=(individus[i]->valeursTraitees[j]-moyennes[i])*
	    (individus[i]->valeursTraitees[j]-moyennes[i]);
	}
      ecartsType[i]=sqrt(ecartsType[i]/(double)(nbDimensions-1));
    }

  for(i=0;i<nbIndividus;i++)
    {
      correlations[i][i]=1.0;
      for(j=i+1;j<nbIndividus;j++)
	{
	  correlations[i][j]=0;
	  for(k=0;k<nbDimensions;k++)
	    {
	      correlations[i][j]+=(individus[i]->valeursTraitees[k]-moyennes[i])*
		(individus[j]->valeursTraitees[k]-moyennes[j]);
	    }
	  correlations[i][j]/=(double)nbDimensions;
	  correlations[j][i]=correlations[i][j];
	}
    }

  for(i=0;i<nbIndividus;i++)
    {
      individus[i]->valeursTraitees=
	(double *)ckrealloc(individus[i]->valeursTraitees,sizeof(double)*nbIndividus);
      for(j=0;j<nbIndividus;j++)
	{
	  individus[i]->valeursTraitees[j]=correlations[i][j];
	  if(individus[i]->valeursTraitees[j]<0)
	    {
	      individus[i]->valeursTraitees[j]=0.0;
	    } 	  
	}
    }

  /*allocation memoire*/
  for(i=0;i<nbIndividus;i++)
    {
      ckfree(correlations[i]);
    }
  ckfree(correlations);
  ckfree(ecartsType);
  ckfree(moyennes);
  /*fin allocation memoire*/
  
}

/*********************************************************************/
/*                                                                   */
/*Procedure de reduction en utilisant un noyau                       */
/*on utilise des individus parmi individus1 pour reduire la dimension*/
/*                                                                   */
/*********************************************************************/

void reductionDimensionParNoyau(int nbIndividus1,int nbDimensions,individu_t **individus1,int nbIndividus2,
				individu_t **individus2,int choixNoyau,int nbPivots,int *newNbDimensions) {
  int i,j;
  /*la taille de la base*/
  int m,*selection; 
  double **alpha1,**alpha2,**Kmm,**Knm1,**Knm2,**covariance,*moyenne;

  m=nbPivots;
  printf("nbPivots : %d %d\n",m,nbIndividus1);

  /*allocation memoire*/
  alpha1=(double **)ckalloc(sizeof(double *)*nbIndividus1);
  alpha2=(double **)ckalloc(sizeof(double *)*nbIndividus2);
  Knm1=(double **)ckalloc(sizeof(double *)*nbIndividus1);
  Knm2=(double **)ckalloc(sizeof(double *)*nbIndividus2);
  Kmm=(double **)ckalloc(sizeof(double *)*m);
  covariance=(double **)ckalloc(sizeof(double *)*m);
  for(i=0;i<nbIndividus1;i++) {
    alpha1[i]=(double *)ckalloc(sizeof(double)*m);
    Knm1[i]=(double *)ckalloc(sizeof(double)*m);
  }
  for(i=0;i<nbIndividus2;i++) {
    alpha2[i]=(double *)ckalloc(sizeof(double)*m);
    Knm2[i]=(double *)ckalloc(sizeof(double)*m);
  }
  for(i=0;i<m;i++) {
    Kmm[i]=(double *)ckalloc(sizeof(double)*m);
    covariance[i]=(double *)ckalloc(sizeof(double)*m);
  }
  selection=(int *)ckalloc(sizeof(int)*m);
  moyenne=(double *)ckalloc(sizeof(double)*m);
   /*fin allocation memoire*/

  /*selection aleatoire de m individus*/
  int nouvelIndividu;
  for(i=0;i<m;i++) {
    while(1) {
      selection[i]=(int)floor((double)nbIndividus1*runif(0,1));
      nouvelIndividu=1;
      for(j=0;j<i-1;j++) {
	if(selection[i]==selection[j]) {
	  nouvelIndividu=0;
	  break;
	}
      }
      if(nouvelIndividu==1) {
	break;
      }
    }
  } 
  
  /*calcul de la matrice Kmm*/
  for(i=0;i<m;i++) {
    for(j=i;j<m;j++) {
      Kmm[i][j]=kernel(nbDimensions,individus1[selection[i]],
		       individus1[selection[j]],choixNoyau);
      Kmm[j][i]=Kmm[i][j];
    }
  }

  /*calcul de la matrice Knm pour 1 et 2*/
  for(i=0;i<nbIndividus1;i++) {
    for(j=0;j<m;j++) {
      Knm1[i][j]=kernel(nbDimensions,individus1[i],individus1[selection[j]],choixNoyau);
    }
  }  

  for(i=0;i<nbIndividus2;i++) {
    for(j=0;j<m;j++) {
      Knm2[i][j]=kernel(nbDimensions,individus2[i],individus1[selection[j]],choixNoyau);
    }
  }

  /*calcul de alpha*/
  inversionMatrice(m,Kmm);

  int k;
  for(i=0;i<nbIndividus1;i++) {
    for(j=0;j<m;j++) {
      alpha1[i][j]=0;
      for(k=0;k<m;k++) {
	alpha1[i][j]+=Knm1[i][k]*Kmm[k][j];
      }
    }
  }
  
  for(i=0;i<nbIndividus2;i++) {
    for(j=0;j<m;j++) {
      alpha2[i][j]=0;
      for(k=0;k<m;k++) {
	alpha2[i][j]+=Knm2[i][k]*Kmm[k][j];
      }
    }
  }
  
  /*calcul de la matrice de variance covariance*/
  calculMatriceCovariance(m,nbIndividus1,alpha1,covariance);

  /*calcul des valeurs propres*/
  double **vecteursPropres,*valeursPropres,*valeursPropres2;
  /*allocation memoire*/
  vecteursPropres=(double **)ckalloc(sizeof(double *)*m);
  for(i=0;i<m;i++) {
    vecteursPropres[i]=(double *)ckalloc(sizeof(double)*m);
  }
  valeursPropres=(double *)ckalloc(sizeof(double)*m);
  valeursPropres2=(double *)ckalloc(sizeof(double)*m);
  /*fin allocation memoire*/
  
  symeig(covariance,m,vecteursPropres,valeursPropres);

  double sommeTemp,sommeValeursPropres=0;
  for(i=0;i<m;i++) {
    printf("%f\n",valeursPropres[i]);
    sommeValeursPropres+=valeursPropres[i];
    valeursPropres2[m-1-i]=valeursPropres[i];
  }
  sommeTemp=0;
  for(i=0;i<m;i++) {
    sommeTemp+=valeursPropres[i];
    if(sommeTemp/sommeValeursPropres>0.9) {
      break;
    }
  }
  *newNbDimensions=i;
  
  clusteringDissimilaritesHolder(m,valeursPropres2,newNbDimensions,0);
  *newNbDimensions=8;
  printf("nouvelle dimension : %d\n",*newNbDimensions);
 

  double **nouvellesCoordonnees1,**nouvellesCoordonnees2;

  /*allocation memoire*/
  nouvellesCoordonnees1=(double **)ckalloc(sizeof(double *)*nbIndividus1);
  nouvellesCoordonnees2=(double **)ckalloc(sizeof(double *)*nbIndividus2);
  for(i=0;i<nbIndividus1;i++) {
    nouvellesCoordonnees1[i]=(double *)ckalloc(sizeof(double)*(*newNbDimensions));
  }
  for(i=0;i<nbIndividus2;i++) {
    nouvellesCoordonnees2[i]=(double *)ckalloc(sizeof(double)*(*newNbDimensions));
  }
  /*fin allocation memoire*/

  /*calcul des nouvelles coordonnees pour les deux ensembles d'individus 1 et 2*/
  for(i=0;i<nbIndividus1;i++) {
    for(j=0;j<*newNbDimensions;j++) {
       nouvellesCoordonnees1[i][j]=0;
       for(k=0;k<m;k++) {
	 nouvellesCoordonnees1[i][j]+=vecteursPropres[k][m-1-j]*(alpha1[i][k]-moyenne[k]);
       }
       nouvellesCoordonnees1[i][j]+=moyenne[j];
    }
  }
  for(i=0;i<nbIndividus2;i++) {
    for(j=0;j<*newNbDimensions;j++) {
       nouvellesCoordonnees2[i][j]=0;
       for(k=0;k<m;k++) {
	 nouvellesCoordonnees2[i][j]+=vecteursPropres[k][m-1-j]*(alpha2[i][k]-moyenne[k]);
       }
       nouvellesCoordonnees2[i][j]+=moyenne[j];
    }
  }
  
  /*mise a jour des coordonnees des individus 1 et 2*/
  for(i=0;i<nbIndividus1;i++) {
    ckfree(individus1[i]->valeursTraitees);
    individus1[i]->valeursTraitees=(double *)ckalloc(sizeof(double)*(*newNbDimensions));
    for(j=0;j<*newNbDimensions;j++) {
      individus1[i]->valeursTraitees[j]=nouvellesCoordonnees1[i][j];
    }
  }
  for(i=0;i<nbIndividus2;i++) {
    ckfree(individus2[i]->valeursTraitees);
    individus2[i]->valeursTraitees=(double *)ckalloc(sizeof(double)*(*newNbDimensions));
    for(j=0;j<*newNbDimensions;j++) {
      individus2[i]->valeursTraitees[j]=nouvellesCoordonnees2[i][j];
    }
  }
  
  
  /*desallocation memoire*/
  for(i=0;i<nbIndividus1;i++) {
    ckfree(alpha1[i]);
    ckfree(Knm1[i]);
  }
  for(i=0;i<nbIndividus2;i++) {
    ckfree(alpha2[i]);
    ckfree(Knm2[i]);
  }
  for(i=0;i<m;i++) {
    ckfree(Kmm[i]);
    ckfree(covariance[i]);
  }
  ckfree(alpha1);
  ckfree(alpha2);
  ckfree(Knm1);
  ckfree(Knm2);
  ckfree(Kmm);
  ckfree(covariance);
  ckfree(selection);
  ckfree(valeursPropres);
  ckfree(valeursPropres2);
  for(i=0;i<nbIndividus1;i++) {
    ckfree(nouvellesCoordonnees1[i]);
  }
  ckfree(nouvellesCoordonnees1);
  for(i=0;i<nbIndividus2;i++) {
    ckfree(nouvellesCoordonnees2[i]);
  }
  ckfree(nouvellesCoordonnees2);
  for(i=0;i<m;i++) {
    ckfree(vecteursPropres[i]);
  }
  ckfree(vecteursPropres);
  /*fin desallocation memoire*/
}
