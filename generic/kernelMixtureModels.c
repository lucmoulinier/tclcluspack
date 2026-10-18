#include "mainClusPack.h"
#include "kernelKmeans.h"
#include "kernelMixtureModels.h"
#include "emKernels.h"
#include "kernels.h"
#include "kmeans.h"
#include "tcl.h"

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/

void kernelMixtureModels(int nbIndividus,
			 int nbDimensions,
			 individu_t **individus,
			 mixtureModel_t *mixtureModel,
			 int nbIterationsMax,
			 int nbClustersSelection,
			 int nbClustersMax,
			 int choixNoyau,
			 int choixDistribution,
			 int accelerationCalculs,
			 int nbPivots,
			 int methodeEstimationParametres,
			 double **kernels) {
  /*declaration de variables*/
  int i,j,k,meilleurNombreClusters,classificationReussie;
  double meilleureVraisemblance,vraisemblance;
  cluster_t *clusters,*meilleursClusters;
  /*fin declaration de variables*/

  if(choixDistribution==EXPONENTIAL) {
    methodeEstimationParametres=MOMENTS;
  }
  
  //printf("debut de kernel mixturemodel\n");
  if(nbClustersMax==-1) {
    nbClustersMax=(int)sqrt((double)nbIndividus);
    if((nbClustersMax<nbIndividus/(nbDimensions+1))&&(nbIndividus/(nbDimensions+1)<60))
      {
	nbClustersMax=nbIndividus/(nbDimensions+1);
      }
  }
  //printf("nbClustersMax : %d\n",nbClustersMax);

  /*allocation memoire*/
  clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*nbClustersMax); 
  meilleursClusters=(cluster_t *)ckalloc(sizeof(cluster_t)*nbClustersMax);
  for(i=0;i<nbClustersMax;i++)
    {
      clusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
      if(choixDistribution==EXPONENTIAL) {
	clusters[i].parametres1D=(double *)ckalloc(sizeof(double));
	meilleursClusters[i].parametres1D=(double *)ckalloc(sizeof(double));
      } else if(choixDistribution==GAMMA) {
	clusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
	meilleursClusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
      } else if(choixDistribution==BETA) {
	clusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
	meilleursClusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
      } else {
	printf("bad choice for distribution\n");
	exit(0);
      }
      if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	clusters[i].alphas=(double *)ckalloc(sizeof(double)*nbIndividus);
	meilleursClusters[i].alphas=(double *)ckalloc(sizeof(double)*nbIndividus);
      }
    }
  for(i=0;i<nbClustersMax;i++)
    {
      meilleursClusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
    }
  for(i=0;i<nbIndividus;i++)
    {
      individus[i]->vrais=(double *)ckalloc(sizeof(double)*nbClustersMax);
    }
  /*fin allocation memoire*/

  if(nbClustersSelection==FIXED)
    {
      /*classification a nombre de clusters fixe*/
      classificationReussie=
	kKernelMixtureModels(nbIndividus,
			     nbDimensions,
			     individus,
			     mixtureModel->nbClusters,
			     clusters,
			     nbIterationsMax,
			     &vraisemblance,
			     choixNoyau,
			     choixDistribution,
			     accelerationCalculs,
			     nbPivots,
			     methodeEstimationParametres,
			     kernels); 
    }
  else
    {
      meilleureVraisemblance=0;
      meilleurNombreClusters=-1;
      
      /*classification*/
      for(i=1;i<nbClustersMax;i++)
	{
	  //printf("cluster teste : %d\n",i);
	  /*classification a nombre de clusters fixe*/
	  classificationReussie=
	    kKernelMixtureModels(nbIndividus,
				 nbDimensions,
				 individus,
				 i,
				 clusters,
				 nbIterationsMax,
				 &vraisemblance,
				 choixNoyau,
				 choixDistribution,
				 accelerationCalculs,
				 nbPivots,
				 methodeEstimationParametres,
				 kernels);
	 
	  if(classificationReussie==OUI)
	    { 
	      if(choixDistribution==EXPONENTIAL) {
		if(nbClustersSelection==AIC) {
		  vraisemblance-=(double)(i-1+i*nbDimensions+i);
		} else {
		  vraisemblance-=0.5*log((double)nbIndividus)*(double)(i-1+i*nbDimensions+i);
		}
	      } else if(choixDistribution==GAMMA) {
		if(nbClustersSelection==AIC) {
		  vraisemblance-=(double)(i-1+i*nbDimensions+2*i);
		} else {
		  vraisemblance-=0.5*log((double)nbIndividus)*(double)(i-1+i*nbDimensions+2*i);
		}
	      } else if(choixDistribution==BETA) {
		if(nbClustersSelection==AIC) {
		  vraisemblance-=(double)(i-1+i*nbDimensions+2*i);
		} else {
		  vraisemblance-=0.5*log((double)nbIndividus)*(double)(i-1+i*nbDimensions+2*i);
		}
	      }
	      //printf("number of clusters : %d ; penalized likelihood : %f\n",i,vraisemblance);

	      if((meilleurNombreClusters==-1)||(vraisemblance>meilleureVraisemblance))
		{
		  meilleureVraisemblance=vraisemblance;
		  meilleurNombreClusters=i;
		  int somme=0;
		  for(j=0;j<meilleurNombreClusters;j++)
		    {
		      meilleursClusters[j].nbIndividus=clusters[j].nbIndividus;
		      if(j<meilleurNombreClusters) {
			somme+=meilleursClusters[j].nbIndividus;
		      }
		      for(k=0;k<clusters[j].nbIndividus;k++)
			{
			  meilleursClusters[j].individus[k]=clusters[j].individus[k];
			}
		    }
		  for(j=0;j<meilleurNombreClusters;j++) {
		    meilleursClusters[j].tau=(double)meilleursClusters[j].nbIndividus/(double)somme;
		  }
		  for(j=0;j<meilleurNombreClusters;j++)
		    {
		      if(choixDistribution==EXPONENTIAL) {
			meilleursClusters[j].parametres1D[0]=clusters[j].parametres1D[0];
			meilleursClusters[j].constanteNormalisation=clusters[j].constanteNormalisation;
		      } else if(choixDistribution==GAMMA) {
			meilleursClusters[j].parametres1D[0]=clusters[j].parametres1D[0];
			meilleursClusters[j].parametres1D[1]=clusters[j].parametres1D[1];
			meilleursClusters[j].constanteNormalisation=clusters[j].constanteNormalisation;
		      } else if(choixDistribution==BETA) {
			meilleursClusters[j].parametres1D[0]=clusters[j].parametres1D[0];
			meilleursClusters[j].parametres1D[1]=clusters[j].parametres1D[1];
		      }
		      meilleursClusters[j].sommeKxixj=clusters[j].sommeKxixj;
		      if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
			for(k=0;k<clusters[j].nbIndividus;k++) {
			  meilleursClusters[j].alphas[k]=clusters[j].alphas[k];
			}
		      }
		    }
		}
	    }
	}
      for(i=0;i<meilleurNombreClusters;i++)
	{
	  for(j=0;j<meilleursClusters[i].nbIndividus;j++)
	    {
	      meilleursClusters[i].individus[j]->cluster=i;
	    }
	}
      //printf("meilleur nombre : %d\n",meilleurNombreClusters);

      mixtureModel->nbClusters=meilleurNombreClusters;

      /*allocation memoire*/
      mixtureModel->clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*mixtureModel->nbClusters);
      for(i=0;i<mixtureModel->nbClusters;i++) {
	if(choixDistribution==EXPONENTIAL) {
	  mixtureModel->clusters[i].parametres1D=(double *)ckalloc(sizeof(double));
	} else if(choixDistribution==GAMMA) {
	  mixtureModel->clusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
	} else if(choixDistribution==BETA) {
	  mixtureModel->clusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
	}
	mixtureModel->clusters[i].qualite=NORMAL;
	mixtureModel->clusters[i].nbIndividus=meilleursClusters[i].nbIndividus;
	mixtureModel->clusters[i].individus=
	  (individu_t **)ckalloc(sizeof(individu_t *)*mixtureModel->clusters[i].nbIndividus);
	mixtureModel->clusters[i].tau=meilleursClusters[i].tau;
	if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	  mixtureModel->clusters[i].alphas=(double *)ckalloc(sizeof(double)*mixtureModel->clusters[i].nbIndividus);
	}
      }
      /*fin allocation memoire*/
      
      for(i=0;i<mixtureModel->nbClusters;i++) {
	for(j=0;j<mixtureModel->clusters[i].nbIndividus;j++) {
	  mixtureModel->clusters[i].individus[j]=meilleursClusters[i].individus[j];
	}
	if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	  for(j=0;j<mixtureModel->clusters[i].nbIndividus;j++) {
	    mixtureModel->clusters[i].alphas[j]=meilleursClusters[i].alphas[j];
	  }
	}
	if(choixDistribution==EXPONENTIAL) {
	  mixtureModel->clusters[i].parametres1D[0]=meilleursClusters[i].parametres1D[0];
	  mixtureModel->clusters[i].constanteNormalisation=meilleursClusters[i].constanteNormalisation;
	} else if(choixDistribution==GAMMA) {
	  mixtureModel->clusters[i].parametres1D[0]=meilleursClusters[i].parametres1D[0];
	  mixtureModel->clusters[i].parametres1D[1]=meilleursClusters[i].parametres1D[1];
	  mixtureModel->clusters[i].constanteNormalisation=meilleursClusters[i].constanteNormalisation;
	} else if(choixDistribution==BETA) {
	  mixtureModel->clusters[i].parametres1D[0]=meilleursClusters[i].parametres1D[0];
	  mixtureModel->clusters[i].parametres1D[1]=meilleursClusters[i].parametres1D[1];
	}
	mixtureModel->clusters[i].sommeKxixj=meilleursClusters[i].sommeKxixj;
      }
    }
  
  /*desallocation memoire*/
  for(i=0;i<nbClustersMax;i++)
    {
      ckfree(clusters[i].individus);
      ckfree(clusters[i].parametres1D);
    }
  ckfree(clusters);
  for(i=0;i<nbClustersMax;i++)
    {
      ckfree(meilleursClusters[i].individus);
    }
  ckfree(meilleursClusters);
  for(i=0;i<nbIndividus;i++)
    {
      ckfree((char *)individus[i]->vrais);
    }  
  /*fin desallocation memoire*/
}


/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*par modeles de melanges                                    */
/*                                                           */
/*************************************************************/

int kKernelMixtureModels(int nbIndividus,
			 int nbDimensions,
			 individu_t **individus,
			 int nbClusters,
			 cluster_t *clusters,
			 int nbIterationsMax,
			 double *meilleureVraisemblance,
			 int choixNoyau,
			 int choixDistribution,
			 int accelerationCalculs,
			 int nbPivots,
			 int methodeEstimationParametres,
			 double **kernels) {
  /*declaration de variables*/
  int i,j,k,meilleurCluster,classificationReussie=NON,compteurIterations=0,resultat,cluster;
  double **points,vraisemblancePrecedente,vraisemblanceCourante,logVraiMax;
  double vraisemblanceIndividu;
  cluster_t *meilleursClusters;
  /*fin declaration de variables*/

  *meilleureVraisemblance=-DBL_MAX;

  //printf("debut kKernel Mixture Models\n");

  /*allocation memoire*/
  points=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      points[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  meilleursClusters=(cluster_t *)ckalloc(sizeof(cluster_t)*nbClusters);
  for(i=0;i<nbClusters;i++)
    {
      meilleursClusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
      if(choixDistribution==EXPONENTIAL) {
	meilleursClusters[i].parametres1D=(double *)ckalloc(sizeof(double));
      } else if(choixDistribution==GAMMA) {
	meilleursClusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
      } else if(choixDistribution==BETA) {
	meilleursClusters[i].parametres1D=(double *)ckalloc(sizeof(double)*2);
      }
      if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	meilleursClusters[i].alphas=(double *)ckalloc(sizeof(double)*nbIndividus);
      }
    }
  /*fin allocation memoire*/

  /*classification initiale par kmeans*/
  resultat=kmeans(nbIndividus,nbDimensions,individus,nbClusters);
  /*
    resultat=kernelKmeans(nbIndividus,
    nbDimensions,individus,nbClusters,choixNoyau,
    accelerationCalculs,kernels);
  */
  //printf("resultat : %d\n",resultat);

  vraisemblancePrecedente=0;
  while(1)
    {
      if(resultat==ECHEC)
	{
	  break;
	}

      for(i=0;i<nbClusters;i++)
	{
	  clusters[i].nbIndividus=0;
	}
      for(i=0;i<nbIndividus;i++)
	{
	  cluster=individus[i]->cluster;
	  clusters[cluster].individus[clusters[cluster].nbIndividus]=individus[i];
	  (clusters[cluster].nbIndividus)++;
	}
      
      for(i=0;i<nbClusters;i++)
	{
	  clusters[i].tau=(double)clusters[i].nbIndividus/(double)nbIndividus;
	  if(clusters[i].nbIndividus<nbDimensions+1)
	    { 
	      compteurIterations=nbIterationsMax+1;
	      break;
	    }
	}
      if(compteurIterations>nbIterationsMax)
	{
	  break;
	}

      /*on maximise la vraisemblance pour chaque cluster*/
      for(i=0;i<nbClusters;i++)
	{
	  for(j=0;j<clusters[i].nbIndividus;j++)
	    {
	      for(k=0;k<nbDimensions;k++)
		{
		  points[j][k]=clusters[i].individus[j]->valeursTraitees[k];
		}
	    }

	  /*etape M : maximum de vraisemblance pour une loi*/
	  if(methodeEstimationParametres==MOMENTS) {
	    estimationParametresKernelMoments(nbDimensions,
					      &(clusters[i]),
					      accelerationCalculs,
					      choixNoyau,
					      choixDistribution,
					      kernels);
	  } else {
	    estimationParametresKernelMaximumLikelihood(nbDimensions,
							&(clusters[i]),accelerationCalculs,choixNoyau,
						       choixDistribution,kernels,-1.0);
	  }
	}

      /*ATTENTION*/
      /*double thetaMoyen=0;
      for(i=0;i<nbClusters;i++) {
	thetaMoyen+=clusters[i].parametres1D[0]*(double)clusters[i].nbIndividus;
      }
      thetaMoyen/=(double)nbIndividus;
      for(i=0;i<nbClusters;i++) {
	clusters[i].parametres1D[0]=thetaMoyen;
	}*/
      /*ATTENTION*/

      /*on calcule la vraisemblance des points par rapport a chaque cluster*/
      calculVraisemblancesKernels(nbIndividus,nbDimensions,individus,nbClusters,clusters,choixNoyau,
				  choixDistribution,methodeEstimationParametres,accelerationCalculs,kernels);
      //printf("Fin Etape M\n");
      
      /*calcul vraisemblance du modele*/
      vraisemblanceCourante=0;
      double maxLogVraisemblance;
      for(i=0;i<nbIndividus;i++)
	{
	  vraisemblanceIndividu=0;
	  maxLogVraisemblance=individus[i]->vrais[0];
	  for(j=0;j<nbClusters;j++)
	    {
	      vraisemblanceIndividu+=exp(individus[i]->vrais[j]);
	      if(individus[i]->vrais[j]>maxLogVraisemblance) {
		maxLogVraisemblance=individus[i]->vrais[j];
	      }
	    } 
	 	  
	  if(!isfinite(log(vraisemblanceIndividu))) {
	    vraisemblanceCourante+=maxLogVraisemblance;
	  } else {
	    vraisemblanceCourante+=log(vraisemblanceIndividu);
	  }
	  //	  printf("vraise courante : %f\n",vraisemblanceCourante);
	}
      
      if((vraisemblanceCourante>*meilleureVraisemblance)||(compteurIterations==0))
	{
	  *meilleureVraisemblance=vraisemblanceCourante;
	  int somme=0;
	  for(i=0;i<nbClusters;i++)
	    { 
	      meilleursClusters[i].nbIndividus=clusters[i].nbIndividus;
	      somme+=meilleursClusters[i].nbIndividus;
	      for(j=0;j<clusters[i].nbIndividus;j++)
		{
		  meilleursClusters[i].individus[j]=clusters[i].individus[j];
		}
	    }
	  for(i=0;i<nbClusters;i++)
	    { 
	      meilleursClusters[i].tau=
		(double)meilleursClusters[i].nbIndividus/(double)somme;
	    } 
	  for(i=0;i<nbClusters;i++)
	    {
	      if(choixDistribution==EXPONENTIAL) {
		meilleursClusters[i].parametres1D[0]=clusters[i].parametres1D[0];
		meilleursClusters[i].constanteNormalisation=clusters[i].constanteNormalisation;
	      } else if(choixDistribution==GAMMA) {
		meilleursClusters[i].parametres1D[0]=clusters[i].parametres1D[0];
		meilleursClusters[i].parametres1D[1]=clusters[i].parametres1D[1];
		meilleursClusters[i].constanteNormalisation=clusters[i].constanteNormalisation;
	      } else if(choixDistribution==BETA) {
		meilleursClusters[i].parametres1D[0]=clusters[i].parametres1D[0];
		meilleursClusters[i].parametres1D[1]=clusters[i].parametres1D[1];
	      }
	      meilleursClusters[i].sommeKxixj=clusters[i].sommeKxixj;
	      if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
		for(j=0;j<clusters[i].nbIndividus;j++) {
		  meilleursClusters[i].alphas[j]=clusters[i].alphas[j];
		}
	      }
	    }
	}

      /*printf("###########################################################\n");
      printf("vraisemblanceCourante : %f (%d)\n",vraisemblanceCourante,compteurIterations);
      printf("###########################################################\n");*/
      
      
      if(fabs(vraisemblanceCourante-vraisemblancePrecedente)<1.0)
	{
	  break;
	}
      vraisemblancePrecedente=vraisemblanceCourante;

      for(i=0;i<nbIndividus;i++)
	{
	  logVraiMax=individus[i]->vrais[0];
	  meilleurCluster=0;
	  for(j=1;j<nbClusters;j++)
	    {
	      if(logVraiMax<individus[i]->vrais[j])
		{
		  logVraiMax=individus[i]->vrais[j];
		  meilleurCluster=j;
		}
	    }
	  individus[i]->cluster=meilleurCluster;
	}
      classificationReussie=OUI;
      compteurIterations++;
    }

  if(classificationReussie==OUI)
    {
      for(i=0;i<nbClusters;i++)
	{ 
	  clusters[i].nbIndividus=meilleursClusters[i].nbIndividus;
	  for(j=0;j<meilleursClusters[i].nbIndividus;j++)
	    {
	      clusters[i].individus[j]=meilleursClusters[i].individus[j];
	    }  
	  if(choixDistribution==EXPONENTIAL) {
	    clusters[i].parametres1D[0]=meilleursClusters[i].parametres1D[0];
	    clusters[i].constanteNormalisation=meilleursClusters[i].constanteNormalisation;	    
	  } else if(choixDistribution==GAMMA) {
	    clusters[i].parametres1D[0]=meilleursClusters[i].parametres1D[0];
	    clusters[i].parametres1D[1]=meilleursClusters[i].parametres1D[1];
	    clusters[i].constanteNormalisation=meilleursClusters[i].constanteNormalisation;	  
	  } else if(choixDistribution==BETA) {
	    clusters[i].parametres1D[0]=meilleursClusters[i].parametres1D[0];
	    clusters[i].parametres1D[1]=meilleursClusters[i].parametres1D[1];
	  }
	  if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	    for(j=0;j<meilleursClusters[i].nbIndividus;j++) {
	      clusters[i].alphas[j]=meilleursClusters[i].alphas[j];
	    }
	  }
	  clusters[i].sommeKxixj=meilleursClusters[i].sommeKxixj;
	}
    }
  
  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++)
    {
      ckfree(points[i]);
    }
  ckfree(points);
  for(i=0;i<nbClusters;i++)
    {
      ckfree(meilleursClusters[i].individus);
      ckfree(meilleursClusters[i].parametres1D);
      if(methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	ckfree(meilleursClusters[i].alphas);
      }
    }
  ckfree(meilleursClusters);
  /*fin desallocation memoire*/;
  
  return classificationReussie;
}

