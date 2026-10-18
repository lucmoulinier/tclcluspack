#include "mainClusPack.h"
#include "mixtureModels.h"
#include "matrice.h"
#include "emGaussian.h"
#include "tcl.h"

/***********************************************************/
/*                                                         */
/*Estimation de la densite par modeles de melanges gaussien*/
/*                                                         */
/***********************************************************/

void mixtureModelDensityEstimation(int nbIndividus,int nbDimensions,individu_t **individus,int diag,
				   mixtureModel_t *mixtureModel) {
  int i;
  int nbIterationsMax=60;

  /*modelisation de l'ensemble des points par un modele de melanges gaussien*/
  int nbClustersMax=(int)sqrt((double)nbIndividus);
  if((nbClustersMax<nbIndividus/(nbDimensions+1))&&(nbIndividus/(nbDimensions+1)<60))
    {
      nbClustersMax=nbIndividus/(nbDimensions+1);
    }

  /*le modele est utilise ensuite : clusteringMixtureModels le conserve
    (ainsi que individus[i]->vrais, alloue par lui) quand destroy vaut 0*/
  mixtureModel->destroy=0;
  clusteringMixtureModels(nbIndividus,nbDimensions,individus,mixtureModel,nbIterationsMax,
			  AIC,GAUSSIAN,diag,nbClustersMax,NON);
  
  printf("nb clusters trouves : %d\n",mixtureModel->nbClusters);
  
  int indicePoint,cluster,meilleurCluster=0;
  double meilleureVraisemblance;
  
  /*affectation des points aux differents clusters du modele de melange*/
  calculVraisemblancesPointsGaussienne(nbIndividus,nbDimensions,individus,
				       mixtureModel->nbClusters-1,mixtureModel->clusters,diag);
  
  for(indicePoint=0;indicePoint<nbIndividus;indicePoint++) { 
    meilleureVraisemblance=-DBL_MAX;
    individus[indicePoint]->vraisemblance=0;
    /*mixtureModel->clusters contient nbClusters-1 elements (le dernier est le bruit)*/
    for(cluster=0;cluster<mixtureModel->nbClusters-1;cluster++) {
      if(mixtureModel->clusters[cluster].qualite==NORMAL) {
	if(individus[indicePoint]->vrais[cluster]>meilleureVraisemblance) {
	  meilleureVraisemblance=individus[indicePoint]->vrais[cluster];
	  meilleurCluster=cluster;
	}
	individus[indicePoint]->vraisemblance+=exp(individus[indicePoint]->vrais[cluster]);
      } 
    }
    individus[indicePoint]->marque=meilleurCluster; 
  }
      

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    ckfree((char *)individus[i]->vrais);
  }
  /*fin desallocation memoire*/
} 

