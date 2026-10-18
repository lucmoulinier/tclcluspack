#include "mainClusPack.h"
#include "kmeans.h"
#include "GAPkmeans.h"

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/
/*le parametre diag indique si l'eventuelle matrice de variance covariance est diagonale*/

void GAPkmeans(int nbIndividus,int nbDimensions,individu_t **individus,int nbClustersMax,
	       int *meilleurNombreClusters) {
  int i,j,k,nbClusters,*clustersDesIndividus;
  double gapCourant,meilleurGap;
  double *valeursMinimum,*valeursMaximum;
  individu_t **individusAleatoires;

  /*desallocation memoire*/
  valeursMinimum=(double *)ckalloc(sizeof(double)*nbDimensions);
  valeursMaximum=(double *)ckalloc(sizeof(double)*nbDimensions);
  individusAleatoires=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      individusAleatoires[i]=(individu_t *)ckalloc(sizeof(individu_t));
      individusAleatoires[i]->id=i;
      individusAleatoires[i]->cluster=AUCUN_GROUPE;
      individusAleatoires[i]->valeursTraitees =(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  clustersDesIndividus=(int *)ckalloc(sizeof(int)*nbIndividus);
  /*fin allocation memoire*/

  if(nbClustersMax==-1) {
    nbClustersMax=(int)sqrt((double)nbIndividus);
    if((nbClustersMax<nbIndividus/(nbDimensions+1))&&(nbIndividus/(nbDimensions+1)<60))
      {
	nbClustersMax=nbIndividus/(nbDimensions+1);
      }
  }

  for(i=0;i<nbDimensions;i++)
    {
      valeursMinimum[i]=individus[0]->valeursTraitees[i];
      valeursMaximum[i]=individus[0]->valeursTraitees[i];
    }
  
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  if(valeursMinimum[j]>individus[i]->valeursTraitees[j])
	    {
	      valeursMinimum[j]=individus[i]->valeursTraitees[j];
	    }
	  if(valeursMaximum[j]<individus[i]->valeursTraitees[j])
	    {
	      valeursMaximum[j]=individus[i]->valeursTraitees[j];
	    }
	}
    }

  int nbSimulations=10,nbSucces,resultatKmeans;
  double logInertieIntraClasseReelle,logInertieIntraClasseAttendu;
  meilleurGap=-DBL_MAX;
  for(nbClusters=1;nbClusters<nbClustersMax;nbClusters++) {
    /*clustering*/
    resultatKmeans=kmeans(nbIndividus,nbDimensions,individus,nbClusters);
    
    if(resultatKmeans==SUCCES) {
      logInertieIntraClasseReelle=
	log(calculInertieIntraClasse(nbIndividus,nbDimensions,individus,nbClusters));
      
      logInertieIntraClasseAttendu=0;
      nbSucces=0;
      for(i=0;i<nbSimulations;i++) {
	/*generation d'un ensemble de points aleatoires dans le pave contenant les points reels*/
	for(j=0;j<nbIndividus;j++) {
	  for(k=0;k<nbDimensions;k++) {
	    individusAleatoires[j]->valeursTraitees[k]=valeursMinimum[k]+
	      runif(0,1)*(valeursMaximum[k]-valeursMinimum[k]);
	  }
	}
	
	/*clustering*/
	resultatKmeans=kmeans(nbIndividus,nbDimensions,individusAleatoires,nbClusters);
	if(resultatKmeans==SUCCES) {
	  logInertieIntraClasseAttendu+=
	    log(calculInertieIntraClasse(nbIndividus,nbDimensions,individusAleatoires,nbClusters));
	  nbSucces++;
	}
      }
      logInertieIntraClasseAttendu/=(double)nbSucces;
      gapCourant=logInertieIntraClasseAttendu-logInertieIntraClasseReelle;
      
      if(gapCourant>meilleurGap) {
	*meilleurNombreClusters=nbClusters;
	meilleurGap=gapCourant;
	for(i=0;i<nbIndividus;i++) {
	  clustersDesIndividus[i]=individus[i]->cluster;
	}
      }
    }
  }
  printf("meilleurNombreDeClusters : %d\n",*meilleurNombreClusters);

  for(i=0;i<nbIndividus;i++) {
    individus[i]->cluster=clustersDesIndividus[i];
  }

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    ckfree(individusAleatoires[i]->valeursTraitees);
    ckfree(individusAleatoires[i]);
  }
  ckfree(individusAleatoires);
  ckfree(valeursMinimum);
  ckfree(valeursMaximum);
  ckfree(clustersDesIndividus);
  /*fin desallocation memoire*/
}


/***************************************************************/
/*                                                             */
/*Procedure de calcul de l'inertie intra-classe d'un clustering*/
/*                                                             */
/***************************************************************/

double calculInertieIntraClasse(int nbIndividus,int nbDimensions,individu_t **individus,
				int nbClusters) {
  int i,j,*nbIndividusParClasse;
  double **centres;

  /*allocation memoire*/
  centres=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    centres[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  nbIndividusParClasse=(int *)ckalloc(sizeof(double)*nbClusters);
  /*fin allocation memoire*/

  for(i=0;i<nbClusters;i++)
    {
      nbIndividusParClasse[i]=0;
    }

  for(i=0;i<nbIndividus;i++)
    {
      nbIndividusParClasse[individus[i]->cluster]++;
    }

  for(i=0;i<nbClusters;i++)
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
  
  for(i=0;i<nbClusters;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  centres[i][j]/=(double)nbIndividusParClasse[i];
	}
    }

  int cluster;
  double temp,inertieIntraClasse=0;
  for(i=0;i<nbIndividus;i++) {
    cluster=individus[i]->cluster;
    for(j=0;j<nbDimensions;j++) {
      temp=individus[i]->valeursTraitees[j]-centres[cluster][j];
      inertieIntraClasse+=temp*temp;
    }
  }

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    ckfree(centres[i]);
  }
  ckfree(centres);
  ckfree(nbIndividusParClasse);
  /*fin desallocation memoire*/

  return inertieIntraClasse;
}
