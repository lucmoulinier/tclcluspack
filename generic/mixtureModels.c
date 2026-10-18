#include "mainClusPack.h"
#include "kmeans.h"
#include "mixtureModels.h"
#include "emGaussian.h"
#include "emSphericalGaussian.h"
#include "emSimplex.h"
#include "tcl.h"

/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*                                                           */
/*************************************************************/
/*le parametre diag indique si l'eventuelle matrice de variance covariance est diagonale*/

int clusteringMixtureModels(int nbIndividus,
			    int nbDimensions,
			    individu_t **individus,
			    mixtureModel_t *mixtureModel,
			    int nbIterationsMax,
			    int nbClustersSelection,
			    int typeDonnees,
			    int diag,
			    int nbClustersMax,
			    int strictExecution)
{
  /*declaration de variables*/
  int i,j,k,l,meilleurNombreClusters,classificationReussie,destroy;
  double logDensiteLoiUniforme,meilleureVraisemblance,vraisemblance;
  cluster_t *clusters,*meilleursClusters;
  /*fin declaration de variables*/
 
  if(diag==OUI) {
    //printf("parameter matrix is diagonal\n");
  }
  //printf("debut de mixturemodel\n");
  if(nbClustersMax==-1) {
    nbClustersMax=(int)sqrt((double)nbIndividus);
    if((nbClustersMax<nbIndividus/(nbDimensions+1))&&(nbIndividus/(nbDimensions+1)<60))
      {
	nbClustersMax=nbIndividus/(nbDimensions+1);
      }
  }
  if(nbClustersSelection==FIXED) {
    nbClustersMax=mixtureModel->nbClusters;
  }
  destroy=mixtureModel->destroy;

  /* allocation memoire*/
  clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*(nbClustersMax+1)); 
  meilleursClusters=(cluster_t *)ckalloc(sizeof(cluster_t)*(nbClustersMax+1));
  for(i=0;i<nbClustersMax+1;i++) {
    clusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
    clusters[i].mu=(double *)ckalloc(sizeof(double)*nbDimensions);
    meilleursClusters[i].mu=(double *)ckalloc(sizeof(double)*nbDimensions);
    
    if((diag==OUI)&&(typeDonnees==GAUSSIAN)) {
      meilleursClusters[i].sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
      clusters[i].sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
    } else {
      clusters[i].sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
      meilleursClusters[i].sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
      for(j=0;j<nbDimensions;j++) {
	clusters[i].sigma[j]=(double *)ckalloc(sizeof(double)*nbDimensions);
	meilleursClusters[i].sigma[j]=(double *)ckalloc(sizeof(double)*nbDimensions);
      }
      if(typeDonnees==SIMPLEX) {
	clusters[i].alphas=(double *)ckalloc(sizeof(double)*nbDimensions);
      }
    }
  }
  for(i=0;i<nbClustersMax+1;i++) {
    meilleursClusters[i].nbIndividus=0;
    meilleursClusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  }
  for(i=0;i<nbIndividus;i++) {
    individus[i]->vrais=(double *)ckalloc(sizeof(double)*(nbClustersMax+1));
  }
  /*fin allocation memoire*/
  
  /*calcul densite uniforme*/
  if((typeDonnees==COORDONNEES)||(typeDonnees==GAUSSIAN)) {
    logDensiteLoiUniforme=calculLogDensiteUniformeGaussienne(nbIndividus,nbDimensions,individus);
  } else if(typeDonnees==SPHERICAL_GAUSSIAN) {
    logDensiteLoiUniforme=calculLogDensiteUniformeGaussienneSpherique(nbIndividus,nbDimensions,individus);
  } else if(typeDonnees==SIMPLEX) {
    logDensiteLoiUniforme=calculLogDensiteUniformeSimplex(nbIndividus,nbDimensions,individus);
  } else {
    /*lm rajoute. ne sert qu'a enleve erreur compilation */
    logDensiteLoiUniforme=0.0;
  }
  
  if(nbClustersSelection==FIXED) {
    /*classification a nombre de clusters fixe*/
    classificationReussie=kmixtureModels(nbIndividus,
					 nbDimensions,
					 individus,
					 mixtureModel->nbClusters,
					 clusters,
					 nbIterationsMax,
					 &vraisemblance,
					 logDensiteLoiUniforme,
					 typeDonnees,
					 diag,
					 strictExecution);
    if(classificationReussie==ECHEC) {
      //fprintf(stderr,"clustering with fixed number of clusters failed\n");
      //fflush(stderr);
      return ECHEC;
    }
    
    for(i=0;i<mixtureModel->nbClusters+1;i++) {
      for(j=0;j<clusters[i].nbIndividus;j++) {
	clusters[i].individus[j]->cluster=i;
      }
    }

    /*clusters 0..nbClusters-1 plus le groupe de bruit d'indice nbClusters :
      la valeur renvoyee (destroy) est meilleurNombreClusters+1*/
    meilleurNombreClusters=mixtureModel->nbClusters;
    if(destroy) {
      for(i=0;i<nbIndividus;i++)
	ckfree((char *)individus[i]->vrais);
      goto done;
    }

    /*allocation memoire*/
    mixtureModel->clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*(mixtureModel->nbClusters+1));
    for(i=0;i<mixtureModel->nbClusters;i++) {
      mixtureModel->clusters[i].mu=(double *)ckalloc(sizeof(double)*nbDimensions);
      if(diag==NON) {
	mixtureModel->clusters[i].sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
	for(j=0;j<nbDimensions;j++) {
	  mixtureModel->clusters[i].sigma[j]=(double *)ckalloc(sizeof(double)*nbDimensions);
	}
      } else {
	mixtureModel->clusters[i].sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
      }
      mixtureModel->clusters[i].qualite=NORMAL;
    }
    for(i=0;i<mixtureModel->nbClusters+1;i++) {
      mixtureModel->clusters[i].nbIndividus=clusters[i].nbIndividus;
      mixtureModel->clusters[i].individus=
	(individu_t **)ckalloc(sizeof(individu_t *)*mixtureModel->clusters[i].nbIndividus);
      for(j=0;j<mixtureModel->clusters[i].nbIndividus;j++) {
	mixtureModel->clusters[i].individus[j]=clusters[i].individus[j];
      }
    }
    /*fin allocation memoire*/
    
    /*le groupe de bruit*/
    if(clusters[mixtureModel->nbClusters].nbIndividus==0) {
      mixtureModel->clusters[mixtureModel->nbClusters].qualite=VIDE;
    } else {
      mixtureModel->clusters[mixtureModel->nbClusters].qualite=BRUIT;
    }
    mixtureModel->clusters[mixtureModel->nbClusters].logDensiteLoiUniforme=logDensiteLoiUniforme;
    
    if(typeDonnees==GAUSSIAN) {
      for(i=0;i<mixtureModel->nbClusters;i++) {
	if(diag==NON) {
	  for(j=0;j<nbDimensions;j++) {
	    mixtureModel->clusters[i].mu[j]=clusters[i].mu[j];
	    for(k=0;k<nbDimensions;k++) {
	      mixtureModel->clusters[i].sigma[j][k]=clusters[i].sigma[j][k];
	    }
	  }
	} else {
	  for(j=0;j<nbDimensions;j++) {
	    mixtureModel->clusters[i].mu[j]=clusters[i].mu[j];
	    mixtureModel->clusters[i].sigmaDiag[j]=clusters[i].sigmaDiag[j];
	  }
	}
	mixtureModel->clusters[i].tau=clusters[i].tau;
	//printf("tau : %f\n",clusters[i].tau);
	if(mixtureModel->clusters[i].tau<ZERO_LIMIT) {
	  mixtureModel->clusters[i].qualite=VIDE;
	}
      }	
    }
  } else {
    /* nombre de cluster a definir */
    meilleureVraisemblance=0;
    meilleurNombreClusters=-1;
    mixtureModel->nbClusters=0;
    
    /* classification */
    for(i=1;i<nbClustersMax;i++) {
      //printf("i : %d\n",i);
      /*classification a nombre de clusters fixe*/
      classificationReussie=kmixtureModels(nbIndividus,
					   nbDimensions,
					   individus,
					   i,
					   clusters,
					   nbIterationsMax,
					   &vraisemblance,
					   logDensiteLoiUniforme,
					   typeDonnees,
					   diag,
					   strictExecution);
      if(classificationReussie==OUI) {
	if((typeDonnees==COORDONNEES)||(typeDonnees==GAUSSIAN)) {
	  /*parametres par cluster : moyenne + covariance (pleine ou diagonale)*/
	  double nbParametres;
	  if(diag==OUI) {
	    nbParametres=(double)(2*nbDimensions)*(double)i+(double)(i+1);
	  } else {
	    nbParametres=(double)((nbDimensions+1)*nbDimensions/2+nbDimensions)*(double)i+
	      (double)(i+1);
	  }
	  /*on maximise logL-penalite : AIC=-2logL+2p, BIC=-2logL+p.ln(n)*/
	  if(nbClustersSelection==AIC) {
	    vraisemblance-=nbParametres;
	  } else {
	    vraisemblance-=0.5*nbParametres*log((double)nbIndividus);
	  }
	} else if(typeDonnees==SPHERICAL_GAUSSIAN) {
	  if(nbClustersSelection==AIC) {
	    vraisemblance-=((double)(nbDimensions-1)*(double)nbDimensions*0.5+
			    (double)nbDimensions-1.0)*(double)i+ (double)(i);
	  } else {
	    vraisemblance-=0.5*(log((double)nbIndividus))*
	      ((double)((nbDimensions-1)*nbDimensions*0.5+nbDimensions-1)*(double)i+ (double)(i) );
	  }
	} else if(typeDonnees==SIMPLEX) {
	  if(nbClustersSelection==AIC) {
	    vraisemblance-=(double)(nbDimensions*i+i);
	  } else {
	    vraisemblance-=(double)(nbDimensions*i+i)*log((double)nbIndividus)*0.5;
	  }
	}
	//printf("number of clusters : %d ; penalized likelihood : %f\n",i,vraisemblance);
	
	if((meilleurNombreClusters==-1)||(vraisemblance>meilleureVraisemblance)) {
	  //printf("meilleur nombre courant : %d\n",meilleurNombreClusters);
	  meilleureVraisemblance=vraisemblance;
	  meilleurNombreClusters=i;
	  int somme=0;
	  for(j=0;j<meilleurNombreClusters+1;j++) {
	    meilleursClusters[j].nbIndividus=clusters[j].nbIndividus;
	    if(j<meilleurNombreClusters) {
	      somme+=meilleursClusters[j].nbIndividus;
	    }
	    if((typeDonnees==GAUSSIAN)&&(j<meilleurNombreClusters)) {
	      if(diag==NON) {
		for(k=0;k<nbDimensions;k++) {
		  meilleursClusters[j].mu[k]=clusters[j].mu[k];
		  for(l=0;l<nbDimensions;l++) {
		    meilleursClusters[j].sigma[k][l]=clusters[j].sigma[k][l];
		  }
		}
	      } else {
		for(k=0;k<nbDimensions;k++) {
		  meilleursClusters[j].mu[k]=clusters[j].mu[k];
		  meilleursClusters[j].sigmaDiag[k]=clusters[j].sigmaDiag[k];
		}
	      }
	    }
	    for(k=0;k<clusters[j].nbIndividus;k++) {
	      meilleursClusters[j].individus[k]=clusters[j].individus[k];
	    }
	  }
	  for(j=0;j<meilleurNombreClusters;j++) {
	    meilleursClusters[j].tau=(double)meilleursClusters[j].nbIndividus/(double)somme;
	  }
	}
      } else {
	/* echec classification */
	//fprintf(stderr,"classification echec\n");
	//fflush(stderr);
      }
    }
    /*fin de boucle sur le nombre de clusters*/

    for(i=0;i<meilleurNombreClusters+1;i++) {
      for(j=0;j<meilleursClusters[i].nbIndividus;j++) {
	meilleursClusters[i].individus[j]->cluster=i;
      }
    }
    
    if(destroy) {
      //ckfree((char *)mixtureModel);
      for(i=0;i<nbIndividus;i++)
	ckfree((char *)individus[i]->vrais);
      goto done;
    }

    /*allocation memoire*/
    mixtureModel->clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*meilleurNombreClusters);
    for(i=0;i<meilleurNombreClusters;i++) {
      mixtureModel->clusters[i].mu=(double *)ckalloc(sizeof(double)*nbDimensions);
      if(diag==NON) {
	mixtureModel->clusters[i].sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
	for(j=0;j<nbDimensions;j++) {
	  mixtureModel->clusters[i].sigma[j]=(double *)ckalloc(sizeof(double)*nbDimensions);
	}
      } else {
	mixtureModel->clusters[i].sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
      }
      mixtureModel->clusters[i].qualite=NORMAL;
    }
    mixtureModel->clusters[meilleurNombreClusters-1].qualite=BRUIT;
    mixtureModel->clusters[meilleurNombreClusters-1].logDensiteLoiUniforme=logDensiteLoiUniforme;
    for(i=0;i<meilleurNombreClusters;i++) {
      mixtureModel->clusters[i].nbIndividus=meilleursClusters[i].nbIndividus;
      mixtureModel->clusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*meilleursClusters[i].nbIndividus);
    }
    /*fin allocation memoire*/
    mixtureModel->nbClusters=meilleurNombreClusters+1;
    
    if(typeDonnees==GAUSSIAN) {
      for(i=0;i<meilleurNombreClusters;i++) {
	for(j=0;j<mixtureModel->clusters[i].nbIndividus;j++) {
	  mixtureModel->clusters[i].individus[j]=meilleursClusters[i].individus[j];
	}
      }
      for(i=0;i<meilleurNombreClusters-1;i++) {
	if(diag==NON) {
	  for(j=0;j<nbDimensions;j++) {
	    mixtureModel->clusters[i].mu[j]=meilleursClusters[i].mu[j];
	    for(k=0;k<nbDimensions;k++) {
	      mixtureModel->clusters[i].sigma[j][k]=meilleursClusters[i].sigma[j][k];
	    }
	  }
	} else {
	  for(j=0;j<nbDimensions;j++) {
	    mixtureModel->clusters[i].mu[j]=meilleursClusters[i].mu[j];
	    mixtureModel->clusters[i].sigmaDiag[j]=meilleursClusters[i].sigmaDiag[j];
	  }
	}
	mixtureModel->clusters[i].tau=meilleursClusters[i].tau;
	if(mixtureModel->clusters[i].tau<ZERO_LIMIT) {
	  mixtureModel->clusters[i].qualite=VIDE;
	}
      }
    }
  }

 done:

  /*desallocation memoire*/
  for(i=0;i<nbClustersMax+1;i++) {
    ckfree((char *)clusters[i].individus);
    ckfree((char *)clusters[i].mu);
    ckfree((char *)meilleursClusters[i].mu);
    /*lm remplace = par == */
    if((diag==OUI)&&(typeDonnees==GAUSSIAN)) {
      ckfree((char *)clusters[i].sigmaDiag); 
      ckfree((char *)meilleursClusters[i].sigmaDiag);
    } else {
      for(j=0;j<nbDimensions;j++) {
	ckfree((char *)clusters[i].sigma[j]); 
	ckfree((char *)meilleursClusters[i].sigma[j]);
      } 
      ckfree((char *)clusters[i].sigma);
      ckfree((char *)meilleursClusters[i].sigma);
    }
    if(typeDonnees==SIMPLEX) {
      ckfree((char *)clusters[i].alphas);
    }
  }
  ckfree((char *)clusters);
  for(i=0;i<nbClustersMax+1;i++) {
    for(j=0;j<nbIndividus;j++) {
      //ckfree((char *)meilleursClusters[i].individus[j]);
    }
    ckfree((char *)meilleursClusters[i].individus);
  }
  ckfree((char *)meilleursClusters);
  /*fin desallocation memoire*/

  if(destroy)
    return meilleurNombreClusters+1;

  return SUCCES;
}


/*************************************************************/
/*                                                           */
/*Procedure de classification en un nombre de classes donnees*/
/*par modeles de melanges                                    */
/*                                                           */
/*************************************************************/

int kmixtureModels(int nbIndividus,
		   int nbDimensions,
		   individu_t **individus,
		   int nbClusters,
		   cluster_t *clusters,
		   int nbIterationsMax,
		   double *meilleureVraisemblance,
		   double logDensiteLoiUniforme,
		   int typeDonnees,
		   int diag,
		   int strictExecution)
{
  /*declaration de variables*/
  int i,j,k,meilleurCluster,classificationReussie=NON,compteurIterations=0,resultat,nbEssais,nbEssaisMax;
  double **points,vraisemblancePrecedente,vraisemblanceCourante,logVraiMax;
  double vraisemblanceIndividu,logVraisemblanceBruit;
  cluster_t *meilleursClusters;
  /*fin declaration de variables*/

  *meilleureVraisemblance=-DBL_MAX;
    
  nbEssaisMax=10;

  /*allocation memoire*/
  points=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    points[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  meilleursClusters=(cluster_t *)ckalloc(sizeof(cluster_t)*(nbClusters+1));
  for(i=0;i<nbClusters+1;i++) {
    meilleursClusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
    meilleursClusters[i].mu=(double *)ckalloc(sizeof(double)*nbDimensions);
    if (diag==NON) {
      meilleursClusters[i].sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
      for(j=0;j<nbDimensions;j++) {
	meilleursClusters[i].sigma[j]=(double *)ckalloc(sizeof(double)*nbDimensions);
      }
    } else {
      meilleursClusters[i].sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
    if (typeDonnees==SIMPLEX) {
      meilleursClusters[i].alphas=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  }
  /*fin allocation memoire*/

  /*classification initiale par kmeans*/
  nbEssais=0;
  while(nbEssais<nbEssaisMax) {
    resultat=kmeans(nbIndividus,nbDimensions,individus,nbClusters);
    if(resultat==SUCCES) {
      for(i=0;i<nbClusters+1;i++)
	{
	  clusters[i].nbIndividus=0;
	}
      for(i=0;i<nbIndividus;i++)
	{
	  clusters[individus[i]->cluster].individus
	    [clusters[individus[i]->cluster].nbIndividus]=individus[i];
	  (clusters[individus[i]->cluster].nbIndividus)++;
	  
	}
      
      for(i=0;i<nbClusters;i++) {
	if((diag==NON)&&(clusters[i].nbIndividus<nbDimensions+1)) { 
	  resultat=ECHEC;
	  break;
	}
      }
    }
    if((strictExecution==NON)||(resultat!=ECHEC)) {
      break;
    }
    nbEssais++;
  } 

  if (resultat==ECHEC) {
    /*desallocation memoire*/
    for(i=0;i<nbIndividus;i++) {
      ckfree((char *)points[i]);
    }
    ckfree((char *)points);
    for(i=0;i<nbClusters+1;i++) {
      ckfree((char *)meilleursClusters[i].individus);
      ckfree((char *)meilleursClusters[i].mu);
      if(diag==NON) {
	for(j=0;j<nbDimensions;j++) {
	  ckfree((char *)meilleursClusters[i].sigma[j]);
	}
	ckfree((char *)meilleursClusters[i].sigma);
      } else {
	ckfree((char *)meilleursClusters[i].sigmaDiag);
      }
      if (typeDonnees==SIMPLEX) {
	ckfree((char *)meilleursClusters[i].alphas);
      }
    }
    ckfree((char *)meilleursClusters);
    /*fin desallocation memoire*/;

    return ECHEC;
  } 

  vraisemblancePrecedente=0;
  while(1)
    {
      if(resultat==ECHEC)
	{
	  break;
	}

      for(i=0;i<nbClusters+1;i++)
	{
	  clusters[i].nbIndividus=0;
	}
      for(i=0;i<nbIndividus;i++)
	{
	  clusters[individus[i]->cluster].individus
	    [clusters[individus[i]->cluster].nbIndividus]=individus[i];
	  (clusters[individus[i]->cluster].nbIndividus)++;
	}
      for(i=0;i<nbClusters;i++)
	{
	  if((clusters[i].nbIndividus<nbDimensions+1)&&(diag==NON))
	    { 
	      compteurIterations=nbIterationsMax+1;
	      break;
	    }
	  clusters[i].tau=(double)clusters[i].nbIndividus/(double)nbIndividus;
	}

      if(compteurIterations>nbIterationsMax)
	{
	  classificationReussie=OUI;
	  break;
	}
      /*on maximise la vraisemblance pour chaque cluster*/
      if(clusters[nbClusters].nbIndividus==0)
	{
	  logVraisemblanceBruit=log(0.01)+logDensiteLoiUniforme;
	}
      else
	{
	  logVraisemblanceBruit=log((double)clusters[nbClusters].nbIndividus/(double)nbIndividus)+
	    logDensiteLoiUniforme;
	}      
      
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
	  if((typeDonnees==COORDONNEES)||(typeDonnees==GAUSSIAN)) {
	      maximisationVraisemblanceGaussienne(points,nbDimensions,clusters[i],diag);
	  } else if(typeDonnees==SPHERICAL_GAUSSIAN) {
	    maximisationVraisemblanceGaussienneSpherique(clusters[i].nbIndividus,points,nbDimensions,clusters[i]);
	  } else if(typeDonnees==SIMPLEX) {
	    maximisationVraisemblanceSimplex(points,nbDimensions,clusters[i]);
	  }
	}

      /*on calcule la vraisemblance des points par rapport a chaque cluster*/
      if((typeDonnees==COORDONNEES)||(typeDonnees==GAUSSIAN)) {
	  calculVraisemblancesPointsGaussienne(nbIndividus,nbDimensions,individus,nbClusters,
					       clusters,diag);
      } else if(typeDonnees==SPHERICAL_GAUSSIAN) {
	calculVraisemblancesPointsGaussienneSpherique(nbIndividus,nbDimensions,individus,nbClusters,clusters);
      } else if(typeDonnees==SIMPLEX) {  
	calculVraisemblancesPointsSimplex(nbIndividus,nbDimensions,individus,nbClusters,clusters);
      } 
      
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
	  vraisemblanceIndividu+=exp(logVraisemblanceBruit);
	  if(logVraisemblanceBruit>maxLogVraisemblance) {
	    maxLogVraisemblance=logVraisemblanceBruit;
	  }
	  
	  if(!isfinite(log(vraisemblanceIndividu))) {
	    vraisemblanceCourante+=maxLogVraisemblance;
	  } else {
	    vraisemblanceCourante+=log(vraisemblanceIndividu);
	  }
	  //printf("vraise courante : %f\n",vraisemblanceCourante);
	}
      
      if((vraisemblanceCourante>*meilleureVraisemblance)||(compteurIterations==0))
	{
	  
	  *meilleureVraisemblance=vraisemblanceCourante;
	  int somme=0;
	  for(i=0;i<nbClusters+1;i++)
	    { 
	      meilleursClusters[i].nbIndividus=clusters[i].nbIndividus;
	      somme+=meilleursClusters[i].nbIndividus;
	      /*lm remplace = par == */
	      if(typeDonnees==GAUSSIAN) {
		if(diag==NON) {
		  for(j=0;j<nbDimensions;j++) {
		    meilleursClusters[i].mu[j]=clusters[i].mu[j];
		    for(k=0;k<nbDimensions;k++) {
		      meilleursClusters[i].sigma[j][k]=clusters[i].sigma[j][k];
		    }
		  }
		} else {
		  for(j=0;j<nbDimensions;j++) {
		    meilleursClusters[i].mu[j]=clusters[i].mu[j];
		    meilleursClusters[i].sigmaDiag[j]=clusters[i].sigmaDiag[j];
		  }
		}
	      }
	      for(j=0;j<clusters[i].nbIndividus;j++)
		{
		  meilleursClusters[i].individus[j]=clusters[i].individus[j];
		}
	    }
	  for(i=0;i<nbClusters+1;i++)
	    { 
	      meilleursClusters[i].tau=
		(double)meilleursClusters[i].nbIndividus/(double)somme;
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
	  if(logVraiMax<logVraisemblanceBruit)
	    {
	      /*probablement du bruit*/
	      meilleurCluster=nbClusters;
	    }
	  individus[i]->cluster=meilleurCluster;
	}
      classificationReussie=OUI;
      compteurIterations++;
    }

  if(classificationReussie==OUI) {
      for(i=0;i<nbClusters+1;i++)
	{ 
	  clusters[i].nbIndividus=meilleursClusters[i].nbIndividus;
	  for(j=0;j<meilleursClusters[i].nbIndividus;j++)
	    {
	      clusters[i].individus[j]=meilleursClusters[i].individus[j];
	    } 
	  clusters[i].tau=meilleursClusters[i].tau;
	  if(clusters[i].tau<ZERO_LIMIT) {
	    clusters[i].qualite=VIDE;
	  }
	}
      for(i=0;i<nbClusters;i++) {
	  if(diag==NON) {
	    for(j=0;j<nbDimensions;j++) {
	      clusters[i].mu[j]=meilleursClusters[i].mu[j];
	      for(k=0;k<nbDimensions;k++) {
		clusters[i].sigma[j][k]=meilleursClusters[i].sigma[j][k];
	      }
	    }
	  } else {
	    for(j=0;j<nbDimensions;j++) {
	      clusters[i].mu[j]=meilleursClusters[i].mu[j];
	      clusters[i].sigmaDiag[j]=meilleursClusters[i].sigmaDiag[j];
	    }
	  }
      }
    }

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    ckfree((char *)points[i]);
  }
  ckfree((char *)points);
  for(i=0;i<nbClusters+1;i++) {
    ckfree((char *)meilleursClusters[i].individus);
    ckfree((char *)meilleursClusters[i].mu);
    if(diag==NON) {
      for(j=0;j<nbDimensions;j++) {
	ckfree((char *)meilleursClusters[i].sigma[j]);
      }
      ckfree((char *)meilleursClusters[i].sigma);
    } else {
      ckfree((char *)meilleursClusters[i].sigmaDiag);
    }
    if (typeDonnees==SIMPLEX) {
      ckfree((char *)meilleursClusters[i].alphas);
    }
  }
  ckfree((char *)meilleursClusters);
  /*fin desallocation memoire*/;
  
  return classificationReussie;
}

