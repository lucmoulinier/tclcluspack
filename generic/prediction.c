#include "mainClusPack.h"
#include "filtering.h"

/*****************************************/
/*                                       */
/*Procedure de prediction sur le test set*/
/*                                       */
/*****************************************/

void prediction(int nbDimensions,int nbClusters,cluster_t *clusters,int tailleTestSet,individu_t **testSet,
		mixtureModel_t **mixtureModels,int diag) {
  /*declaration de variables*/
  int i,j,k,meilleurCluster;
  /*fin declaration de variables*/

  printf("PREDICTION\n");
  
  /*calcul des bandwidth pour le point teste courant*/
  int indicePoint,cluster;
  double vraisemblance,meilleureVraisemblance,temp1,temp2;
  
  //filtreDimensions(tailleTestSet,nbDimensions,&nbDimensionsMM,testSet);
  //nbDimensions=nbDimensionsMM;
  
  /*affectation des points aux differents clusters du modele de melange*/
  if(diag==NON) {
    for(indicePoint=0;indicePoint<tailleTestSet;indicePoint++) {
      meilleureVraisemblance=-DBL_MAX;
      for(cluster=0;cluster<nbClusters;cluster++) {
	vraisemblance=0;
	for(k=0;k<mixtureModels[cluster]->nbClusters;k++) {
	  if(mixtureModels[cluster]->clusters[k].qualite==NORMAL) {
	    temp1=0;
	    for(i=0;i<nbDimensions;i++) {
	      temp2=testSet[indicePoint]->valeursTraitees[i]-mixtureModels[cluster]->clusters[k].mu[i];
	      for(j=i+1;j<nbDimensions;j++) {
		temp1+=temp2*(testSet[indicePoint]->valeursTraitees[j]-
			      mixtureModels[cluster]->clusters[k].mu[j])*
		  mixtureModels[cluster]->clusters[k].sigma[i][j];
	      }
	      temp1+=0.5*temp2*temp2*mixtureModels[cluster]->clusters[k].sigma[i][i];
	    }
	    vraisemblance+=clusters[cluster].tau*mixtureModels[cluster]->clusters[k].tau*
	      exp(-temp1-0.5*mixtureModels[cluster]->clusters[k].logDeterminant)/
	      pow(2.0*M_PI,(double)nbDimensions*0.5);
	  }
	}
	if(vraisemblance>meilleureVraisemblance) {
	  meilleureVraisemblance=vraisemblance;
	  meilleurCluster=cluster;
	}
      }
      testSet[indicePoint]->cluster=meilleurCluster;
    }
  } else {
    for(indicePoint=0;indicePoint<tailleTestSet;indicePoint++) {
      meilleureVraisemblance=-DBL_MAX;
      for(cluster=0;cluster<nbClusters;cluster++) {
	vraisemblance=0;
	for(k=0;k<mixtureModels[cluster]->nbClusters;k++) {
	  if(mixtureModels[cluster]->clusters[k].qualite==NORMAL) {
	    temp1=0;
	    for(i=0;i<nbDimensions;i++) {
	      temp2=testSet[indicePoint]->valeursTraitees[i]-mixtureModels[cluster]->clusters[k].mu[i];
	      temp1+=0.5*temp2*temp2*mixtureModels[cluster]->clusters[k].sigmaDiag[i];
	    }
	    vraisemblance+=clusters[cluster].tau*mixtureModels[cluster]->clusters[k].tau*
	      exp(-temp1-0.5*mixtureModels[cluster]->clusters[k].logDeterminant)/
	      pow(2.0*M_PI,(double)nbDimensions*0.5);
	  }
	}
	if(vraisemblance>meilleureVraisemblance) {
	  meilleureVraisemblance=vraisemblance;
	  meilleurCluster=cluster;
	}
      }
      testSet[indicePoint]->cluster=meilleurCluster;
    }
  }
}
