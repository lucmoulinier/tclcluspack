#include "mainClusPack.h"
#include "kernels.h"
#include "emKernels.h"
#include "matrice.h"
#include "filtering.h"
#include "densityEstimation.h"
#include "kNearestNeighbours.h"
#include "symeig.h"
#include "kernelMixtureModels.h"
#include "mixtureModels.h"

/******************************************************************/
/*                                                                */
/*procedure d'apprentissage : k au sens des k plus proches voisins*/
/*                                                                */
/******************************************************************/

void apprentissageKernels(int *nbDimensions,int *tailleTrainingSet,individu_t **trainingSet,int tailleTestSet,
			  individu_t **testSet,int *nbClusters,cluster_t **clusters,mixtureModel_t ***mixtureModels,
			  mixtureModel_t *mixtureModelPourDensite,int bandwidthType,double bandwidth,
			  int accelerationCalculs,int *diag,int choixNoyau,int *choixDistribution,
			  double scale,int methodeEstimationParametres,int verbose,double regularisation) {
  int i,j,cluster;
  int nbDimensionsMax=20;  

  /*recherche du nombre de clusters differents*/
  *nbClusters=-1;
  for(i=0;i<*tailleTrainingSet;i++) {
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

  for(i=0;i<*tailleTrainingSet;i++) {
    (*clusters)[trainingSet[i]->cluster].nbIndividus++;
  }

  /*allocation memoire*/
  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*(*clusters)[i].nbIndividus);
  } 
  *mixtureModels=(mixtureModel_t **)ckalloc(sizeof(mixtureModel_t *)*(*nbClusters));
  for(i=0;i<*nbClusters;i++) {
    (*mixtureModels)[i]=(mixtureModel_t *)ckalloc(sizeof(mixtureModel_t));
  }
  /*fin allocation memoire*/
  
  /*constitution des differents clusters*/
  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].nbIndividus=0;
  }
  for(i=0;i<*tailleTrainingSet;i++) {
    cluster=trainingSet[i]->cluster;
    (*clusters)[cluster].individus[(*clusters)[cluster].nbIndividus]=trainingSet[i];
    (*clusters)[cluster].nbIndividus++;
    /*allocation memoire*/
    trainingSet[i]->inverseBandwidth=(double *)ckalloc(sizeof(double)*(*nbDimensions));
    trainingSet[i]->vrais=(double *)ckalloc(sizeof(double));
    /*fin allocation memoire*/
  }

  if(*nbDimensions>nbDimensionsMax) {
    *diag=OUI;
  } else {
    *diag=NON;
  }
  *diag=OUI;
  /*modelisation de l'ensemble des points par un modele de melanges gaussien*/
   if((bandwidthType==ABMM)||(bandwidthType==ABNN)) { 
     /*estimation de la densite*/    
     
     //filtreDimensionsDeuxEnsembles(*tailleTrainingSet,*nbDimensions,&nbDimensionsMM,trainingSet,tailleTestSet,testSet);
     //*nbDimensions=nbDimensionsMM;
     printf("nbDimension a jour : %d\n",*nbDimensions);
     if(bandwidthType==ABMM) {
 	    /*estimation de la densite par modele de melanges*/
    	    mixtureModelDensityEstimation(*tailleTrainingSet,*nbDimensions,trainingSet,*diag,mixtureModelPourDensite);
    	    printf("fin d'estimation de la densite\n");
  	} 
  }

  for(i=0;i<*nbClusters;i++) {
    (*clusters)[i].tau=(double)(*clusters)[i].nbIndividus/(double)*tailleTrainingSet;
  }

  int kPlusProchesVoisins=10;
  /*calcul du bandwidth*/
  if(bandwidthType==FIXED_BANDWIDTH) {
    calculFixedBandwidth(*nbDimensions,*tailleTrainingSet,trainingSet,bandwidth);
  } else if(bandwidthType==ABMM) {
    /*utilise par defaut*/
    calculMMDensityBandwidth(*nbDimensions,*tailleTrainingSet,trainingSet);
  }  else if(bandwidthType==ABNN) {
    /*utilise par defaut*/
    calculNNAdaptiveBandwidth(*nbDimensions,*tailleTrainingSet,trainingSet,
			      kPlusProchesVoisins,scale,accelerationCalculs,0);
  }  else if(bandwidthType==AB) {
    /*utilise par defaut*/
     calculAdaptiveBandwidth(*nbDimensions,*tailleTrainingSet,trainingSet,kPlusProchesVoisins,scale);
  } else if(bandwidthType==ABCW) {
    /*utilise par defaut*/
     calculAdaptiveBandwidthComponentWise(*nbDimensions,*tailleTrainingSet,trainingSet,kPlusProchesVoisins,scale);
  } 

  /*calcul des lois pour chaque cluster*/
  printf("nbClusters : %d\n",*nbClusters);

  //int nbIterationsMax=60;
  
  double *parametresExponentielles,*parametresGammas,*sommesKxixjSansLambdas,*sommesKxixjAvecLambdas;
  /*allocation memoire*/
  parametresExponentielles=(double *)ckalloc(sizeof(double)*(*nbClusters));
  parametresGammas=(double *)ckalloc(sizeof(double)*(*nbClusters*2));
  sommesKxixjSansLambdas=(double *)ckalloc(sizeof(double)*(*nbClusters));
  sommesKxixjAvecLambdas=(double *)ckalloc(sizeof(double)*(*nbClusters));
   /*fin allocation memoire*/


  for(cluster=0;cluster<*nbClusters;cluster++) {
    printf("cluster %d\n",cluster);
    /*int nbClustersMax;   
    nbClustersMax=(int)sqrt((double)(*clusters)[cluster].nbIndividus);
    if((nbClustersMax<(*clusters)[cluster].nbIndividus/(*nbDimensions+1))&&
       ((*clusters)[cluster].nbIndividus/(*nbDimensions+1)<60)) {
	nbClustersMax=(*clusters)[cluster].nbIndividus/(*nbDimensions+1);
      }*/

    double **kernels;
    if(accelerationCalculs==OUI) {
      kernels=(double **)ckalloc(sizeof(double *)*(*clusters)[cluster].nbIndividus);
      for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
	kernels[i]=(double *)ckalloc(sizeof(double)*(*clusters)[cluster].nbIndividus);
      }
      for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
	(*clusters)[cluster].individus[i]->id=i;
	for(j=i;j<(*clusters)[cluster].nbIndividus;j++) {
	  kernels[i][j]=kernel(*nbDimensions,(*clusters)[cluster].individus[i],
			       (*clusters)[cluster].individus[j],choixNoyau);
	  kernels[j][i]=kernels[i][j];
	}
      }
    }

    (*mixtureModels)[cluster]->nbClusters=1;
    (*mixtureModels)[cluster]->clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*1);
    (*mixtureModels)[cluster]->clusters[0].tau=1.0;
    (*mixtureModels)[cluster]->clusters[0].qualite=NORMAL;
    (*mixtureModels)[cluster]->clusters[0].nbIndividus=(*clusters)[cluster].nbIndividus;
    (*mixtureModels)[cluster]->clusters[0].individus=(individu_t **)ckalloc(sizeof(individu_t *)*(*clusters)[cluster].nbIndividus);
    (*mixtureModels)[cluster]->clusters[0].alphas=(double *)ckalloc(sizeof(double)*(*clusters)[cluster].nbIndividus);
    for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
      (*mixtureModels)[cluster]->clusters[0].individus[i]=(*clusters)[cluster].individus[i];
    }
    if((*clusters)[cluster].nbIndividus>0) {
	    if(regularisation<0) {
	    /*distribution exponentielle*/
	    estimationParametresKernelMoments(*nbDimensions,&((*mixtureModels)[cluster]->clusters[0]),NON,
    			choixNoyau,EXPONENTIAL,kernels);
	    parametresExponentielles[cluster]=(*mixtureModels)[cluster]->clusters[0].parametres1D[0];
	    sommesKxixjSansLambdas[cluster]=(*mixtureModels)[cluster]->clusters[0].sommeKxixj;
	    printf("exp : %f\n",parametresExponentielles[cluster]);
	    }

	    /*distribution gamma*/
	    estimationParametresKernelMaximumLikelihood(*nbDimensions,&((*mixtureModels)[cluster]->clusters[0]),NON,
    			choixNoyau,GAMMA,kernels,regularisation);
	    parametresGammas[2*cluster]=(*mixtureModels)[cluster]->clusters[0].parametres1D[0];
	    parametresGammas[2*cluster+1]=(*mixtureModels)[cluster]->clusters[0].parametres1D[1];
	    sommesKxixjAvecLambdas[cluster]=(*mixtureModels)[cluster]->clusters[0].sommeKxixj;
	   printf("gamma : %f %f\n",parametresGammas[cluster*2],parametresGammas[cluster*2+1]);

    } else {
  	 (*mixtureModels)[cluster]->clusters[0].qualite=VIDE;
    }
   

    if(accelerationCalculs==OUI) {
      for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
	ckfree(kernels[i]);
      }
      ckfree(kernels);
    }
  }
	
 
  if(regularisation>0) {
	*choixDistribution=GAMMA;
  } else {
  /*comparaison des differents modeles par rapport de vraisemblance*/
  int degreExponentielleThetaConstant=1;
  int degreExponentielle=*nbClusters;
  //int degreGamma=*nbClusters*2;

  double logVraisemblanceExponentielleThetaConstant,logVraisemblanceExponentielle,logVraisemblanceGamma;
  
  /*calcul logvraisemblance thetaConstant*/
  double thetaMoyen=0;
  for(cluster=0;cluster<*nbClusters;cluster++) {
  	 thetaMoyen+=parametresExponentielles[cluster]*(double)(*mixtureModels)[cluster]->clusters[0].nbIndividus;
  }
  thetaMoyen/=(double)*tailleTrainingSet;
  logVraisemblanceExponentielleThetaConstant=0;
  for(cluster=0;cluster<*nbClusters;cluster++) {
  	(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=thetaMoyen;
	(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjSansLambdas[cluster];
        for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
        	calculVraisemblancesKernels(1,*nbDimensions,&((*clusters)[cluster].individus[i]),1,
				 &((*mixtureModels)[cluster]->clusters[0]),choixNoyau,EXPONENTIAL,
				 MOMENTS,accelerationCalculs,NULL);
                logVraisemblanceExponentielleThetaConstant+=(*clusters)[cluster].individus[i]->vrais[0];
        }
  } 
  printf("logvraisemblance theta constant : %f\n",logVraisemblanceExponentielleThetaConstant);

  /*calcul logvraisemblance exponentielle*/
  logVraisemblanceExponentielle=0;
  for(cluster=0;cluster<*nbClusters;cluster++) {
  	(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresExponentielles[cluster];
 	(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjSansLambdas[cluster];
        for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
        	calculVraisemblancesKernels(1,*nbDimensions,&((*clusters)[cluster].individus[i]),1,
				 &((*mixtureModels)[cluster]->clusters[0]),choixNoyau,EXPONENTIAL,
				 MOMENTS,accelerationCalculs,NULL);
                logVraisemblanceExponentielle+=(*clusters)[cluster].individus[i]->vrais[0];
        }
  } 
  printf("logvraisemblance exponentielle : %f\n",logVraisemblanceExponentielle);
 
  /*calcul logvraisemblance gamma*/
  logVraisemblanceGamma=0;
  for(cluster=0;cluster<*nbClusters;cluster++) {
  	(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresGammas[2*cluster];
        (*mixtureModels)[cluster]->clusters[0].parametres1D[1]=parametresGammas[2*cluster+1];
	(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjAvecLambdas[cluster];
        for(i=0;i<(*clusters)[cluster].nbIndividus;i++) {
        	calculVraisemblancesKernels(1,*nbDimensions,&((*clusters)[cluster].individus[i]),1,
				 &((*mixtureModels)[cluster]->clusters[0]),choixNoyau,GAMMA,
				 MAXIMUM_LIKELIHOOD,accelerationCalculs,NULL);
                logVraisemblanceGamma+=(*clusters)[cluster].individus[i]->vrais[0];
        }
  } 
  printf("logvraisemblance gamma : %f\n",logVraisemblanceGamma);
 
  /*comparaison theta fixe/exponentielle*/
  double statistique,risqueObserve;
  statistique=2.0*(logVraisemblanceExponentielle-logVraisemblanceExponentielleThetaConstant);
  risqueObserve=pchisq(statistique,degreExponentielle-degreExponentielleThetaConstant,0,0);
  printf("risque observe (test1) : %f\n",risqueObserve);

  FILE *fichierChoixModele;
  
  /*choix par BIC*/
  double BICFixedExponential,BICExponential,BICGamma;
  BICFixedExponential=logVraisemblanceExponentielleThetaConstant-0.5*(double)(*nbClusters*2+1)*log((double)*tailleTrainingSet);
  BICExponential=logVraisemblanceExponentielle-0.5*(double)(*nbClusters*3)*log((double)*tailleTrainingSet);
  BICGamma=logVraisemblanceGamma-0.5*(double)(*nbClusters*4)*log((double)*tailleTrainingSet);

  if((BICFixedExponential>BICExponential)&&(BICFixedExponential>BICGamma)) {
	*choixDistribution=EXPONENTIAL;
	for(cluster=0;cluster<*nbClusters;cluster++) {
		(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=thetaMoyen;
		(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjSansLambdas[cluster];
	}
	if(verbose==OUI) {
		fichierChoixModele=fopen("choixModele.txt","a");
		fprintf(fichierChoixModele,"0\n");
		fclose(fichierChoixModele);	
	}
  } else if ((BICExponential>BICFixedExponential)&&(BICExponential>BICGamma)) {
       	*choixDistribution=EXPONENTIAL;
	for(cluster=0;cluster<*nbClusters;cluster++) {
		(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresExponentielles[cluster];
		(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjSansLambdas[cluster];
	}
	if(verbose==OUI) {
		fichierChoixModele=fopen("choixModele.txt","a");
		fprintf(fichierChoixModele,"1\n");
		fclose(fichierChoixModele);	
	}
  } else {
	*choixDistribution=GAMMA;
	for(cluster=0;cluster<*nbClusters;cluster++) {
		(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresGammas[2*cluster];
		(*mixtureModels)[cluster]->clusters[0].parametres1D[1]=parametresGammas[2*cluster+1];		
		(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjAvecLambdas[cluster];		
	}	
	if(verbose==OUI) {
		fichierChoixModele=fopen("choixModele.txt","a");
		fprintf(fichierChoixModele,"2\n");
		fclose(fichierChoixModele);	
	}	
  }}
 
 /*choix par test loglikelihood ratio 
if(risqueObserve>risqueLimite) {
 	comparaison theta fixe/gamma
	statistique=2.0*(logVraisemblanceGamma-logVraisemblanceExponentielleThetaConstant);
  	risqueObserve=pchisq(statistique,degreGamma-degreExponentielleThetaConstant,0,0);
	if(risqueObserve>risqueLimite) {
		*choixDistribution=EXPONENTIAL;
		for(cluster=0;cluster<*nbClusters;cluster++) {
			(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=thetaMoyen;
			(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjSansLambdas[cluster];
		}
		if(verbose==OUI) {
			fichierChoixModele=fopen("choixModele.txt","a");
			fprintf(fichierChoixModele,"0\n");
			fclose(fichierChoixModele);	
		}
	} else {
		*choixDistribution=GAMMA;
		for(cluster=0;cluster<*nbClusters;cluster++) {
			(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresGammas[2*cluster];
			(*mixtureModels)[cluster]->clusters[0].parametres1D[1]=parametresGammas[2*cluster+1];		
			(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjAvecLambdas[cluster];		
		}	
		if(verbose==OUI) {
			fichierChoixModele=fopen("choixModele.txt","a");
			fprintf(fichierChoixModele,"2\n");
			fclose(fichierChoixModele);	
		}	
	}
  } else {
 	comparaison exponentielle/gamma
 	statistique=2.0*(logVraisemblanceGamma-logVraisemblanceExponentielle);
  	risqueObserve=pchisq(statistique,degreGamma-degreExponentielle,0,0);
	  printf("risque observe (test1) : %f\n",risqueObserve);
   	if(risqueObserve>risqueLimite) {
       		*choixDistribution=EXPONENTIAL;
		for(cluster=0;cluster<*nbClusters;cluster++) {
			(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresExponentielles[cluster];
			(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjSansLambdas[cluster];
		}
		if(verbose==OUI) {
			fichierChoixModele=fopen("choixModele.txt","a");
			fprintf(fichierChoixModele,"1\n");
			fclose(fichierChoixModele);	
		}
  	} else {
		*choixDistribution=GAMMA;
		for(cluster=0;cluster<*nbClusters;cluster++) {
			(*mixtureModels)[cluster]->clusters[0].parametres1D[0]=parametresGammas[2*cluster];
			(*mixtureModels)[cluster]->clusters[0].parametres1D[1]=parametresGammas[2*cluster+1];		
			(*mixtureModels)[cluster]->clusters[0].sommeKxixj=sommesKxixjAvecLambdas[cluster];
		}	
		if(verbose==OUI) {
			fichierChoixModele=fopen("choixModele.txt","a");
			fprintf(fichierChoixModele,"2\n");
			fclose(fichierChoixModele);	
		}
	}
  }}*/
 

  /*desallocation memoire*/
  ckfree(parametresExponentielles);
  ckfree(parametresGammas);
  ckfree(sommesKxixjSansLambdas);
  ckfree(sommesKxixjAvecLambdas);
  for(i=0;i<*tailleTrainingSet;i++) {
    ckfree(trainingSet[i]->vrais);
  }
  /*fin desallocation memoire*/
}

