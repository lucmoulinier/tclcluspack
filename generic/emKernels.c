#include "mainClusPack.h"
#include "kernels.h"
#include "emKernels.h"
#include "emGaussian.h"

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calculVraisemblancesKernels(int nbIndividus,int nbDimensions,individu_t **individus,int nbClusters,
				 cluster_t *clusters,int choixNoyau,int choixDistribution,
				 int methodeEstimationParametres,int accelerationCalculs,double **kernels) {
  int i,j,k;
  double x;

  for(i=0;i<nbClusters;i++)
    {
      /*calcul des log-vraisemblances*/
      for(j=0;j<nbIndividus;j++)
	{
	  individus[j]->vrais[i]=log(clusters[i].tau);
	  if(accelerationCalculs==NON) {
	    if(methodeEstimationParametres==MOMENTS) {
	      x=0;
	      for(k=0;k<clusters[i].nbIndividus;k++) {
		x-=kernel(nbDimensions,individus[j],clusters[i].individus[k],choixNoyau);
	      }
	      
	      x*=2.0/(double)clusters[i].nbIndividus;
	      x+=kernel(nbDimensions,individus[j],individus[j],choixNoyau)+clusters[i].sommeKxixj/
		(double)(clusters[i].nbIndividus*clusters[i].nbIndividus);
	    } else if (methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	      x=0;
	      for(k=0;k<clusters[i].nbIndividus;k++) {
		x-=clusters[i].alphas[k]*kernel(nbDimensions,individus[j],clusters[i].individus[k],choixNoyau);
	      }
	      x*=2.0;
	      x+=kernel(nbDimensions,individus[j],individus[j],choixNoyau)+clusters[i].sommeKxixj;
	    }
	  } else {
	    if(methodeEstimationParametres==MOMENTS) {
	      x=0;
	      for(k=0;k<clusters[i].nbIndividus;k++) {
	        x-=kernels[individus[j]->id][clusters[i].individus[k]->id];
	      }
	    
	      x*=2.0/(double)clusters[i].nbIndividus;
	      x+=kernels[individus[j]->id][individus[j]->id]+clusters[i].sommeKxixj/
	        (double)(clusters[i].nbIndividus*clusters[i].nbIndividus);
	    } else if (methodeEstimationParametres==MAXIMUM_LIKELIHOOD) {
	      x=0;
	      for(k=0;k<clusters[i].nbIndividus;k++) {
		x-=clusters[i].alphas[k]*kernels[individus[j]->id][clusters[i].individus[k]->id];
	      }
	      x*=2.0;
	      x+=kernels[individus[j]->id][individus[j]->id]+clusters[i].sommeKxixj;
	    }
	  }
	  if(choixDistribution==EXPONENTIAL) { 
	    individus[j]->vrais[i]+=-log(clusters[i].parametres1D[0])-x/clusters[i].parametres1D[0];
	    //individus[j]->vrais[i]-=log(clusters[i].constanteNormalisation);
	  } else if(choixDistribution==GAMMA) {
	    individus[j]->vrais[i]+=(clusters[i].parametres1D[1]-1.0)*log(x)-x/clusters[i].parametres1D[0]-
	      clusters[i].parametres1D[1]*log(clusters[i].parametres1D[0])-lgammafn(clusters[i].parametres1D[1]);
	    //individus[j]->vrais[i]-=log(clusters[i].constanteNormalisation);
	  } else if(choixDistribution==BETA) {
	    printf("%f %f\n",clusters[i].parametres1D[0],clusters[i].parametres1D[1]);
	    individus[j]->vrais[i]+=(clusters[i].parametres1D[0]-1.0)*log(x/2.0)+
	      (clusters[i].parametres1D[1]-1.0)*log(1.0-x/2.0)+
	      lgammafn(clusters[i].parametres1D[0]+clusters[i].parametres1D[1])-
	      lgammafn(clusters[i].parametres1D[0])-lgammafn(clusters[i].parametres1D[1]);
	  }
	}
     
    }
}

/***********************************************************************/
/*                                                                     */
/*Procedure de calcul des parametres de la loi (exponentielle ou gamma)*/
/*qui renvoie la log-vraisemblance                                     */
/*                                                                     */
/***********************************************************************/

void estimationParametresKernelMoments(int nbDimensions,cluster_t *cluster,int accelerationCalculs,
				       int choixNoyau,int choixDistribution,double **kernels) {
  int i,j;
  double **points;

  /*allocation memoire*/
  points=(double **)ckalloc(sizeof(double *)*cluster->nbIndividus);
  for(i=0;i<cluster->nbIndividus;i++) {
    points[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  cluster->sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
  cluster->mu=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  for(i=0;i<cluster->nbIndividus;i++) {
    for(j=0;j<nbDimensions;j++) {
      points[i][j]=cluster->individus[i]->valeursTraitees[j];
    }
  }

  /*calcul de la moyenne et de la variance*/
  maximisationVraisemblanceGaussienne(points,nbDimensions,*cluster,OUI);

  /*calcul de la somme des Kxixj*/
  cluster->sommeKxixj=0;
  if(accelerationCalculs==NON) {
    for(i=0;i<cluster->nbIndividus;i++) {
      cluster->sommeKxixj+=kernel(nbDimensions,cluster->individus[i],cluster->individus[i],choixNoyau);
      for(j=i+1;j<cluster->nbIndividus;j++) {
	cluster->sommeKxixj+=2.0*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
      }
    }
  } else {
    for(i=0;i<cluster->nbIndividus;i++) {
      cluster->sommeKxixj+=kernels[cluster->individus[i]->id][cluster->individus[i]->id];
      for(j=i+1;j<cluster->nbIndividus;j++) {
	cluster->sommeKxixj+=2.0*kernels[cluster->individus[i]->id][cluster->individus[j]->id];
      }
    }
  }

  double temp,moyenne,variance,distanceCarre,sommeLog,sommeValeurs;

  /*calcul de la moyenne des distances au carre*/
  moyenne=0;
  for(i=0;i<cluster->nbIndividus;i++) {
    moyenne+=kernel(nbDimensions,cluster->individus[i],cluster->individus[i],choixNoyau);
  }

  moyenne-=cluster->sommeKxixj/(double)cluster->nbIndividus;
  moyenne/=(double)cluster->nbIndividus;

  sommeLog=0;	
  sommeValeurs=0;

  /*calcul de la variance des distances au carre*/
  //  printf("#############################################################\n");
  variance=0;
  if(accelerationCalculs==NON) {
    for(i=0;i<cluster->nbIndividus;i++) {
      distanceCarre=kernel(nbDimensions,cluster->individus[i],cluster->individus[i],choixNoyau);
      temp=0;
      for(j=0;j<cluster->nbIndividus;j++) {
	temp+=kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
      }
      distanceCarre-=2.0*temp/(double)cluster->nbIndividus-cluster->sommeKxixj/
	(double)(cluster->nbIndividus*cluster->nbIndividus);
      //  printf("%f\n",distanceCarre);
      variance+=(distanceCarre-moyenne)*(distanceCarre-moyenne);
      sommeLog+=log(distanceCarre);
      sommeValeurs+=distanceCarre;
    }
    variance/=(double)cluster->nbIndividus;
  } else {
    for(i=0;i<cluster->nbIndividus;i++) {
      distanceCarre=kernels[cluster->individus[i]->id][cluster->individus[i]->id];
      temp=0;
      for(j=0;j<cluster->nbIndividus;j++) {
	temp+=kernels[cluster->individus[i]->id][cluster->individus[j]->id];
      }
      distanceCarre-=2.0*temp/(double)cluster->nbIndividus-cluster->sommeKxixj/
	(double)(cluster->nbIndividus*cluster->nbIndividus);
      variance+=(distanceCarre-moyenne)*(distanceCarre-moyenne);
      sommeLog+=log(distanceCarre);
      sommeValeurs+=distanceCarre;
    }
    variance/=(double)cluster->nbIndividus;
  }

  /*allocation memoire*/
  if(choixDistribution==EXPONENTIAL) {
    cluster->parametres1D=(double *)ckalloc(sizeof(double));
  } else if ((choixDistribution==GAMMA)||(choixDistribution==BETA)) {
    cluster->parametres1D=(double *)ckalloc(sizeof(double)*2);
  } else {
    printf("bad distribution choice\n");
    exit(0);
  }
  /*fin allocation memoire*/
  
  if(choixDistribution==EXPONENTIAL) {
    /*theta*/
    cluster->parametres1D[0]=moyenne;
  } else if(choixDistribution==GAMMA) {
    /*theta*/
    cluster->parametres1D[0]=variance/moyenne;
    /*k*/
    cluster->parametres1D[1]=moyenne*moyenne/variance;
  } else if(choixDistribution==BETA) {
    double theta0,theta1,theta2,m1,m2;
    moyenne/=2.0;
    variance/=4.0;
    m1=moyenne;
    m2=1-moyenne;
    
    /*theta0*/
    theta0=m1*m2/variance-1.0;

    /*theta1*/
    theta1=m1*theta0;
    
    /*theta2*/
    theta2=m2*theta0;
    
    cluster->parametres1D[0]=theta1;
    cluster->parametres1D[1]=theta2;
  } 

  /*on ajuste la constante de normalisation*/
  //cluster->constanteNormalisation=normalisationKernel(nbDimensions,choixNoyau,choixDistribution,cluster);

  /*desallocation memoire*/
  for(i=0;i<cluster->nbIndividus;i++) {
    ckfree(points[i]);
  }
  ckfree(points);
  ckfree(cluster->sigmaDiag);
  ckfree(cluster->mu);
  /*fin desallocation memoire*/
}

/*********************************************************/
/*                                                       */
/*Procedure d'ajustement de la constante de normalisation*/
/*                                                       */
/*********************************************************/

double normalisationKernel(int nbDimensions,int choixNoyau,int choixDistribution,cluster_t *cluster) {
  int i,j,nbIterations;
  individu_t individuAleatoire;
  double x,logRatio;

  /*allocation memoire*/
  individuAleatoire.valeursTraitees=(double *)ckalloc(sizeof(double)*nbDimensions);
  individuAleatoire.inverseBandwidth=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  for(i=0;i<nbDimensions;i++) {
    individuAleatoire.inverseBandwidth[i]=cluster->individus[0]->inverseBandwidth[i];
  }

  nbIterations=10000;
  double constante=0;
  for(i=0;i<nbIterations;i++) {
    for(j=0;j<nbDimensions;j++) {
      individuAleatoire.valeursTraitees[j]=rnorm(cluster->mu[j],sqrt(cluster->sigmaDiag[j]));
    }
    x=0;
    for(j=0;j<cluster->nbIndividus;j++) {
      x-=kernel(nbDimensions,&individuAleatoire,cluster->individus[j],choixNoyau);
    }
    
    x*=2.0/(double)cluster->nbIndividus;
    x+=kernel(nbDimensions,&individuAleatoire,&individuAleatoire,choixNoyau)+cluster->sommeKxixj/
      (double)(cluster->nbIndividus*cluster->nbIndividus);

    if(choixDistribution==EXPONENTIAL) {
      logRatio=-log(cluster->parametres1D[0])-x/cluster->parametres1D[0]; 
    } else {
      logRatio=(cluster->parametres1D[1]-1.0)*log(x)-x/cluster->parametres1D[0]-
	cluster->parametres1D[1]*log(cluster->parametres1D[0])-lgammafn(cluster->parametres1D[1]);
    }

    /*on retranche la log vraisemblance de la gaussienne*/
    for(j=0;j<nbDimensions;j++) {
      logRatio-=dnorm(individuAleatoire.valeursTraitees[j],cluster->mu[j],sqrt(cluster->sigmaDiag[j]),1);
    }
    constante+=exp(logRatio);
  }
  constante/=(double)nbIterations;

  /*desallocation memoire*/
  ckfree(individuAleatoire.valeursTraitees);
  ckfree(individuAleatoire.inverseBandwidth);
  /*fin desallocation memoire*/

  return constante;
}

/***********************************************************************/
/*                                                                     */
/*Procedure de calcul des parametres de la loi (exponentielle ou gamma)*/
/*qui renvoie la log-vraisemblance                                     */
/*                                                                     */
/***********************************************************************/

void estimationParametresKernelMaximumLikelihood(int nbDimensions,cluster_t *cluster,int accelerationCalculs,
					         int choixNoyau,int choixDistribution,double **kernels,
						 double regularisation) {
  int i,j,k;
  double *parameterDerivatives,*alphaDerivatives;

  /*allocation memoire*/
  parameterDerivatives=(double *)ckalloc(sizeof(double)*2);
  alphaDerivatives=(double *)ckalloc(sizeof(double)*cluster->nbIndividus);
  /*fin allocation memoire*/
  
  /*on commence par une estimation par moments*/
  estimationParametresKernelMoments(nbDimensions,cluster,accelerationCalculs,choixNoyau,
				    choixDistribution,kernels);
  
  double temp;
  temp=1.0/(double)cluster->nbIndividus;
  for(i=0;i<cluster->nbIndividus;i++) {
    cluster->alphas[i]=temp;
  }

  cluster->sommeKxixj=0;
  for(i=0;i<cluster->nbIndividus;i++) {
    for(j=0;j<cluster->nbIndividus;j++) {
	cluster->sommeKxixj+=cluster->alphas[i]*cluster->alphas[j]*
	  kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
    }
  }

  int nbDeltas,iteration,nbIterations;
  double delta,*deltas,deltaMax,meilleurDelta,somme,premiereSomme,deuxiemeSomme,*normesCarre,*alphasK;
  double *anciensParametres,*KxixjIfixe,logVraisemblance,meilleureLogVraisemblance,*sauvegardeParametres1D;
  double *sauvegardeAlphas;
  nbIterations=10;
  
  /*allocation memoire*/
  normesCarre=(double *)ckalloc(sizeof(double)*cluster->nbIndividus);
  alphasK=(double *)ckalloc(sizeof(double)*cluster->nbIndividus);
  anciensParametres=(double *)ckalloc(sizeof(double)*2);
  KxixjIfixe=(double *)ckalloc(sizeof(double)*cluster->nbIndividus);
  sauvegardeParametres1D=(double *)ckalloc(sizeof(double)*2);
  sauvegardeAlphas=(double *)ckalloc(sizeof(double)*cluster->nbIndividus);
  /*fin allocation memoire*/
printf("sommeKxixj : %f\n",cluster->sommeKxixj);
  printf("taille : %d\n",cluster->nbIndividus);

  if(choixDistribution==GAMMA) {    
    nbDeltas=10;
    
    /*allocation memoire*/
    deltas=(double *)ckalloc(sizeof(double)*nbDeltas);
    /*fin allocation memoire*/

    deltas[0]=1.0;
    for(i=1;i<nbDeltas;i++) {
	deltas[i]=deltas[i-1]/10.0;
    }
  }

  for(iteration=0;iteration<nbIterations;iteration++) {
    /*********************/
    /*calcul des derivees*/
    /*********************/
    if((choixDistribution==GAMMA)||(choixDistribution==BETA)) {
      for(i=0;i<2;i++) {
	anciensParametres[i]=cluster->parametres1D[i];
      }
      for(i=0;i<cluster->nbIndividus;i++) {
	alphasK[i]=0;
	KxixjIfixe[i]=0;
	for(j=0;j<cluster->nbIndividus;j++) {
	  alphasK[i]+=cluster->alphas[j]*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
	  KxixjIfixe[i]+=kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
	}
	
	normesCarre[i]=0;
	for(j=0;j<cluster->nbIndividus;j++) {
	  normesCarre[i]-=cluster->alphas[j]*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
	}
	normesCarre[i]*=2.0;
	normesCarre[i]+=kernel(nbDimensions,cluster->individus[i],cluster->individus[i],choixNoyau)+cluster->sommeKxixj;
      }
    }
    
    if(choixDistribution==GAMMA) {
      premiereSomme=0;
      deuxiemeSomme=0;
      for(i=0;i<cluster->nbIndividus;i++) {
	//printf("normeCarre %d : %f\n",i,normesCarre[i]);
	premiereSomme+=log(normesCarre[i]);
	deuxiemeSomme+=normesCarre[i];
      }
//      printf("sommes : %f %f\n",premiereSomme,deuxiemeSomme);
      /*theta*/
      parameterDerivatives[0]=deuxiemeSomme/(cluster->parametres1D[0]*cluster->parametres1D[0])-
	cluster->parametres1D[1]*(double)cluster->nbIndividus/cluster->parametres1D[0];
      /*k*/
      parameterDerivatives[1]=premiereSomme-(double)cluster->nbIndividus*
	(log(cluster->parametres1D[0])+digamma(cluster->parametres1D[1]));
      
      /*alphas*/
      for(i=0;i<cluster->nbIndividus;i++) {
	alphaDerivatives[i]=0;
	for(j=0;j<cluster->nbIndividus;j++) {
	  alphaDerivatives[i]+=(-2.0*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau)+
				2.0*alphasK[i])/normesCarre[j];
	}
	alphaDerivatives[i]*=anciensParametres[1]-1;
	alphaDerivatives[i]+=-2.0*(double)cluster->nbIndividus*alphasK[i]/cluster->parametres1D[0]+
	  2.0/cluster->parametres1D[0]*KxixjIfixe[i];

	if(regularisation>0) {
		//alphaDerivatives[i]+=2.0*regularisation*cluster->alphas[i];	
		if(cluster->alphas[i]>0) {
			alphaDerivatives[i]-=regularisation;
		} else {
			alphaDerivatives[i]+=regularisation;
		}
	}
      }
     // printf("parametres deriv : %f %f\n",parameterDerivatives[0],parameterDerivatives[1]);
    } else if(choixDistribution==BETA) {
      premiereSomme=0;
      deuxiemeSomme=0;
      for(i=0;i<cluster->nbIndividus;i++) {
	premiereSomme+=log(normesCarre[i]*0.5);
	deuxiemeSomme+=log(1.0-normesCarre[i]*0.5);
      }
      
      /*theta1*/
      parameterDerivatives[0]=(double)cluster->nbIndividus*digamma(cluster->parametres1D[0]+premiereSomme);
      /*theta2*/
      parameterDerivatives[1]=(double)cluster->nbIndividus*digamma(cluster->parametres1D[1]+deuxiemeSomme);
      //printf("parametres : %f %f\n",parameterDerivatives[0],parameterDerivatives[1]);
      /*alphas*/
      for(i=0;i<cluster->nbIndividus;i++) {
	alphaDerivatives[i]=0;
	premiereSomme=0;
	deuxiemeSomme=0;
	for(j=0;j<cluster->nbIndividus;j++) {
	  premiereSomme+=(-2.0*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau)+
			  2.0*alphasK[i])/normesCarre[j];
	  deuxiemeSomme+=(2.0*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau)-
			  2.0*alphasK[i])/(1-0.5*normesCarre[j]);
	}
	alphaDerivatives[i]=(anciensParametres[0]-1)*premiereSomme+(anciensParametres[1]-1)*deuxiemeSomme;
	//printf("alpha %d : %f\n",i,cluster->alphas[i]);
      }
    }
    
    /****************************/
    /*mise a jour des parametres*/
    /****************************/
    
    /*les alphas et betas doivent rester positifs, cela fixe le delta maximum*/
    deltaMax=DBL_MAX;
    //printf("parametres %f ; %f\n",-cluster->parametres1D[0],cluster->parametres1D[1]);
    for(i=0;i<2;i++) {
      temp=-cluster->parametres1D[i]/parameterDerivatives[i];
      if((temp>0)&&(temp<deltaMax)) {
	//printf("div %f\n",temp);
	deltaMax=temp;
      }
    }
    /*TEMPORAIREfor(i=0;i<cluster->nbIndividus;i++) {
      temp=-cluster->alphas[i]/alphaDerivatives[i];
      if((temp>0)&&(temp<deltaMax)) {
	deltaMax=temp;
      }
      temp=(1.0-cluster->alphas[i])/alphaDerivatives[i];
      if((temp>0)&&(temp<deltaMax)) {
	deltaMax=temp;
      }
    }*/


    meilleureLogVraisemblance=-DBL_MAX;
    /*sauvegarde des parametres courants*/
    for(i=0;i<2;i++) {
       sauvegardeParametres1D[i]=cluster->parametres1D[i];
    }
    for(i=0;i<cluster->nbIndividus-1;i++) {
       sauvegardeAlphas[i]=cluster->alphas[i];
    }
    //printf("deltaMax : %f\n",deltaMax);
//getchar();
    meilleurDelta=deltas[nbDeltas-1];
    for(i=nbDeltas-1;i>=0;i--) {
      delta=deltas[i]; 	
      if(delta>deltaMax) {
      	break;
      }	      

      /*mise a jour des parametres provisoire*/
      for(j=0;j<2;j++) {
        cluster->parametres1D[j]=sauvegardeParametres1D[j]+delta*parameterDerivatives[j];
      }
	/*printf("deriv1 : %f ; deriv2 : %f\n",parameterDerivatives[0],parameterDerivatives[1]);
	printf("Temp : para1 : %f ; para2 : %f\n",cluster->parametres1D[0],cluster->parametres1D[1]);
      */somme=0;
      for(j=0;j<cluster->nbIndividus-1;j++) {
	cluster->alphas[j]=sauvegardeAlphas[j]+delta*alphaDerivatives[j];
        somme+=cluster->alphas[j];
      }
      cluster->alphas[cluster->nbIndividus-1]=1.0-somme;	

      cluster->sommeKxixj=0;
      for(j=0;j<cluster->nbIndividus;j++) {
        for(k=0;k<cluster->nbIndividus;k++) {
	  cluster->sommeKxixj+=cluster->alphas[j]*cluster->alphas[k]*
	    kernel(nbDimensions,cluster->individus[j],cluster->individus[k],choixNoyau);
        }
      }

      logVraisemblance=calculLogVraisemblanceGamma(nbDimensions,cluster,choixNoyau,regularisation);
      //printf("delta : %f ; logvraisemblance : %f\n",delta,logVraisemblance);
      if(logVraisemblance>meilleureLogVraisemblance) {
	meilleureLogVraisemblance=logVraisemblance;
	meilleurDelta=delta;
	printf("meilleur delta : %f\n",meilleurDelta);
      }
    }
    /*mise a jour des parametres definitive pour l'iteration*/
    for(i=0;i<2;i++) {
      cluster->parametres1D[i]=sauvegardeParametres1D[i]+meilleurDelta*parameterDerivatives[i];
    }     
//   printf("para1 : %f ; para2 : %f\n",cluster->parametres1D[0],cluster->parametres1D[1]);
 
    somme=0;
    for(i=0;i<cluster->nbIndividus-1;i++) {
      cluster->alphas[i]=sauvegardeAlphas[i]+meilleurDelta*alphaDerivatives[i];
      somme+=cluster->alphas[i];
	//printf("alpha %d : %f\n",i,cluster->alphas[i]);
    }
    cluster->alphas[cluster->nbIndividus-1]=1.0-somme;
    //printf("alpha %d : %f\n",i,cluster->alphas[i]);
    cluster->sommeKxixj=0;
    for(i=0;i<cluster->nbIndividus;i++) {
      for(j=0;j<cluster->nbIndividus;j++) {
	cluster->sommeKxixj+=cluster->alphas[i]*cluster->alphas[j]*
	  kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
      }
    }
    //printf("sommeKxixj : %f\n",cluster->sommeKxixj);
  }
  /*desallocation memoire*/
  ckfree(parameterDerivatives);
  ckfree(alphaDerivatives);
  ckfree(alphasK);
  ckfree(KxixjIfixe);
  ckfree(normesCarre);
  ckfree(anciensParametres);
  ckfree(sauvegardeParametres1D);
  ckfree(sauvegardeAlphas);
  if(choixDistribution==GAMMA) {
    ckfree(deltas);
  }
  /*fin desallocation memoire*/
}

/**************************************************/
/*                                                */
/* Calcul de la log vraisemblance pour la loi gamma*mei */
/**************************************************/

double calculLogVraisemblanceGamma(int nbDimensions,cluster_t *cluster,int choixNoyau,double regularisation) {
  int i,j;
  double x,logVraisemblance=0;

  /*calcul des log-vraisemblances*/
  for(i=0;i<cluster->nbIndividus;i++) {
    x=0;
    for(j=0;j<cluster->nbIndividus;j++) { 
      x-=cluster->alphas[j]*kernel(nbDimensions,cluster->individus[i],cluster->individus[j],choixNoyau);
    }
    x*=2.0;
    x+=kernel(nbDimensions,cluster->individus[i],cluster->individus[i],choixNoyau)+cluster->sommeKxixj;

    logVraisemblance+=(cluster->parametres1D[1]-1.0)*log(x)-x/cluster->parametres1D[0]-
      cluster->parametres1D[1]*log(cluster->parametres1D[0])-lgammafn(cluster->parametres1D[1]);
  }

  double temp=0;
  if(regularisation>0) {
	for(i=0;i<cluster->nbIndividus;i++) {
		//temp-=cluster->alphas[i]*cluster->alphas[i];
		temp+=fabs(cluster->alphas[i]);
	}
  	logVraisemblance-=temp*regularisation;
  }
  return logVraisemblance;
}

