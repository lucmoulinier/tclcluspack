#include "mainClusPack.h"
#include "apprentissage.h"
#include "matrice.h"
#include "mixtureModels.h"
#include "filtering.h"

/******************************************************************/
/*                                                                */
/*procedure d'apprentissage : k au sens des k plus proches voisins*/
/*                                                                */
/******************************************************************/

void apprentissage(int *nbDimensions,int tailleTrainingSet,individu_t **trainingSet,int tailleTestSet,
		   individu_t **testSet,int *nbClusters,cluster_t **clusters,mixtureModel_t ***mixtureModel,
		   int *diag) {
  int i,j,cluster;
  int nbDimensionsMax=40;

  /*recherche du nombre de clusters differents*/
  *nbClusters=-1;
  for(i=0;i<tailleTrainingSet;i++) {
    if(trainingSet[i]->cluster>*nbClusters) {
      *nbClusters=trainingSet[i]->cluster;
    }
  }
  (*nbClusters)++;
 
  /*allocation memoire*/
  *clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*(*nbClusters));
  /*fin allocation memoire*/

  /*calcul de la taille des differents clusters*/
  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].nbIndividus=0;
  }

  for(i=0;i<tailleTrainingSet;i++) {
    (*clusters)[trainingSet[i]->cluster].nbIndividus++;
  }
  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].tau=(double)(*clusters)[i].nbIndividus/(double)tailleTrainingSet;
  }

  /*allocation memoire*/
  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*(*clusters)[i].nbIndividus);
  }
  *mixtureModel=(mixtureModel_t **)ckalloc(sizeof(mixtureModel_t *)*(*nbClusters));
  for(i=0;i<*nbClusters;i++) {
    (*mixtureModel)[i]=(mixtureModel_t *)ckalloc(sizeof(mixtureModel_t));
  }
  /*fin allocation memoire*/

  /*constitution des differents clusters*/
  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].nbIndividus=0;
  }
  for(i=0;i<tailleTrainingSet;i++) {
    cluster=trainingSet[i]->cluster;
    (*clusters)[cluster].individus[(*clusters)[cluster].nbIndividus]=trainingSet[i];
    (*clusters)[cluster].nbIndividus++;
  }
  
  int nbIterationsMax=60;
  
  /*modelisation de chaque cluster par un modele de melanges gaussien*/
  int nbDimensionsMM;
  double *sigmaDiag;
  filtreDimensionsDeuxEnsembles(tailleTrainingSet,*nbDimensions,&nbDimensionsMM,
				trainingSet,tailleTestSet,testSet);
  *nbDimensions=nbDimensionsMM;
  printf("nbDimensions apres filtre : %d\n",*nbDimensions);
  
  for(cluster=0;cluster<*nbClusters;cluster++) {
    if(*nbDimensions>nbDimensionsMax) {
      *diag=OUI; 
      /*allocation memoire*/
      sigmaDiag=(double *)ckalloc(sizeof(double)*(*nbDimensions));
      /*fin allocation memoire*/
    } else {
      *diag=NON;
    }
    int nbClustersMax;   
    nbClustersMax=(int)sqrt((double)(*clusters)[cluster].nbIndividus);
    if((nbClustersMax<(*clusters)[cluster].nbIndividus/(*nbDimensions+1))&&
       ((*clusters)[cluster].nbIndividus/(*nbDimensions+1)<60))
      {
	nbClustersMax=(*clusters)[cluster].nbIndividus/(*nbDimensions+1);
      }
    (*mixtureModel)[cluster]->destroy=0;
    clusteringMixtureModels((*clusters)[cluster].nbIndividus,*nbDimensions,(*clusters)[cluster].individus,
			    (*mixtureModel)[cluster],nbIterationsMax,AIC,GAUSSIAN,*diag,nbClustersMax,NON);
    //printf("nb clusters trouves : %d\n",(*mixtureModel)[cluster]->nbClusters);
    
    /*inversion de la matrice de variance-covariance et recuperation du determinant*/
    for(i=0;i<(*mixtureModel)[cluster]->nbClusters;i++) {
      if((*mixtureModel)[cluster]->clusters[i].qualite==NORMAL) {
	if(*diag==NON) {
	  (*mixtureModel)[cluster]->clusters[i].logDeterminant=
	    inversionMatrice(*nbDimensions,(*mixtureModel)[cluster]->clusters[i].sigma);
	} else {
	  /*on copie la matrice de covariance dans sigma*/
	  for(j=0;j<*nbDimensions;j++)
	    {
	      sigmaDiag[j]=(*mixtureModel)[cluster]->clusters[i].sigmaDiag[j];
	      printf("sigma %d : %f\n",j,sigmaDiag[j]);
	      if(sigmaDiag[j]==0) {
		sigmaDiag[j]=1.0;
	      }
	    }
	  (*mixtureModel)[cluster]->clusters[i].logDeterminant=0;
	  for(j=0;j<*nbDimensions;j++) {
	    (*mixtureModel)[cluster]->clusters[i].logDeterminant+=log(sigmaDiag[j]);
	    (*mixtureModel)[cluster]->clusters[i].sigmaDiag[j]=1.0/sigmaDiag[j];
	  }
	}
      }
    }
    if(*diag==OUI) {
      ckfree(sigmaDiag);
    }
  }
}
