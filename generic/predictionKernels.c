#include "mainClusPack.h"
#include "kernels.h"
#include "emKernels.h"
#include "emGaussian.h"

/*****************************************/
/*                                       */
/*Procedure de prediction sur le test set*/
/*                                       */
/*****************************************/

void predictionKernels(int nbDimensions,int nbClusters,cluster_t *clusters,int tailleTrainingSet,individu_t **trainingSet,
		       int tailleTestSet,individu_t **testSet,
		       mixtureModel_t **mixtureModels,mixtureModel_t mixtureModelPourDensite,int bandwidthType,
		       double bandwidth,int accelerationCalculs,int diag,int choixNoyau,int choixDistribution,
		       int methodeEstimationParametres) {
  /*declaration de variables*/
  int i,j,k,meilleurCluster;
  /*fin declaration de variables*/

  printf("PREDICTION\n");

  /*allocation memoire*/
  for(i=0;i<tailleTestSet;i++) {
    testSet[i]->vrais=NULL;
    testSet[i]->inverseBandwidth=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  /*fin allocation memoire*/

  if(bandwidthType==FIXED_BANDWIDTH) {
    calculFixedBandwidth(nbDimensions,tailleTestSet,testSet,bandwidth);
  } else if(bandwidthType==ABMM) {
    for(i=0;i<tailleTestSet;i++) {
      testSet[i]->vrais=(double *)ckalloc(sizeof(double)*mixtureModelPourDensite.nbClusters);
    }
    calculVraisemblancesPointsGaussienne(tailleTestSet,nbDimensions,testSet,mixtureModelPourDensite.nbClusters-1,
					 mixtureModelPourDensite.clusters,diag);

    int cluster;
    for(i=0;i<tailleTestSet;i++) { 
      testSet[i]->vraisemblance=0;
      for(cluster=0;cluster<mixtureModelPourDensite.nbClusters;cluster++) {
	if(mixtureModelPourDensite.clusters[cluster].qualite==NORMAL) {
	  testSet[i]->vraisemblance+=exp(testSet[i]->vrais[cluster]);
	} 
      }
    }
    calculMMDensityBandwidth(nbDimensions,tailleTestSet,testSet);
  } else if((bandwidthType==ABNN)||(bandwidthType==AB)) {
    /*pour chaque point, on cherche le point le plus proche pour lui attribuer le meme bandwidth*/
    double distance,distanceMin,temp;
    individu_t *pointLePlusProche;
    for(i=0;i<tailleTestSet;i++) {
	distanceMin=DBL_MAX; 	
	for(j=0;j<tailleTrainingSet;j++) {
    		distance=0;
		for(k=0;k<nbDimensions;k++) {
			temp=testSet[i]->valeursTraitees[k]-trainingSet[j]->valeursTraitees[k];
			distance+=temp*temp;		
		}
		if(distance<distanceMin) {
			distanceMin=distance;	
			pointLePlusProche=trainingSet[j];	
		}   	
	}
	for(k=0;k<nbDimensions;k++) {
	    testSet[i]->inverseBandwidth[k]=pointLePlusProche->inverseBandwidth[k];	
	}
    }
  }
  
  int indicePoint,cluster;
  double vraisemblance,meilleureVraisemblance;
  /*nouvelle version*/
  for(indicePoint=0;indicePoint<tailleTestSet;indicePoint++) {
    meilleureVraisemblance=-DBL_MAX;
    for(cluster=0;cluster<nbClusters;cluster++) {
      /*allocation memoire*/
      if(testSet[indicePoint]->vrais!=NULL) {
	ckfree(testSet[indicePoint]->vrais);
      } 
      testSet[indicePoint]->vrais=(double *)ckalloc(sizeof(double)*mixtureModels[cluster]->nbClusters);
      /*desallocation memoire*/
      
      vraisemblance=0;
      /*calcul des log vraisemblances pour chaque cluster du melange*/

      if(choixDistribution==GAMMA) {
	 methodeEstimationParametres=MAXIMUM_LIKELIHOOD;
      } else {
	 methodeEstimationParametres=MOMENTS;
      }
      calculVraisemblancesKernels(1,nbDimensions,&(testSet[indicePoint]),mixtureModels[cluster]->nbClusters,
	 			     mixtureModels[cluster]->clusters,choixNoyau,choixDistribution,
				     methodeEstimationParametres,0,NULL);
     
      for(k=0;k<mixtureModels[cluster]->nbClusters;k++) {
	//printf("k : %d\n",k);
	if(mixtureModels[cluster]->clusters[k].qualite==NORMAL) {
	  //	  printf("exp : %f\n",exp(testSet[indicePoint]->vrais[k]));
	  vraisemblance+=clusters[cluster].tau*mixtureModels[cluster]->clusters[k].tau*exp(testSet[indicePoint]->vrais[k]);
	}
      }
     
      if(vraisemblance>meilleureVraisemblance) {
	meilleureVraisemblance=vraisemblance;
	meilleurCluster=cluster;
      }
    }
    /*if(meilleureVraisemblance/somme<0.90) {
	printf("%s %f\n",testSet[indicePoint]->nom,vraisemblance/somme);
      }*/
    testSet[indicePoint]->cluster=meilleurCluster;
  }
}









