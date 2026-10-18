#include "mainClusPack.h"
#include "matrice.h"
#include "kernels.h"
#include "space.h"

/**************************************************************/
/*                                                            */
/*Procedure de calcul du noyau par appel aux autres procedures*/
/*                                                            */
/**************************************************************/

double kernel(int nbDimensions,individu_t *point1,individu_t *point2,int choixNoyau) {
  double valeur;
  if(choixNoyau==RADIAL_GAUSSIAN) {
    valeur=radialGaussianKernel(nbDimensions,point1,point2);
  } else if (choixNoyau==EPANECHNIKOV) {
    valeur=EpanechnikovKernel(nbDimensions,point1,point2);
  } 
  else {
    valeur=polynomialKernel(nbDimensions,point1,point2,choixNoyau);
  }
   // printf("ker : %f\n",valeur);
  return valeur;
}

/******************************************************************/
/*                                                                */
/*Procedure de calcul du noyau gaussien                           */
/*bandwidth contient les variances pour les differentes dimensions*/
/*                                                                */
/******************************************************************/
/*noyau special Landsat satellite*/
/*double radialGaussianKernel(int nbDimensions,individu_t *point1,individu_t *point2) {
  double valeur;
  int i,j;
  double temp,valeur1,valeur2;
  double valeurMax;

  valeurMax=0;
  for(j=0;j<4;j++) { 
    valeur1=0;
    valeur2=0;
    for(i=0;i<nbDimensions-4;i++) {
      temp=point1->valeursTraitees[i]-point2->valeursTraitees[(i/4)*4+(i+j)%4];
      temp*=temp;
      valeur1+=temp*point1->inverseBandwidth[i]/2.0;
      valeur2+=temp*point2->inverseBandwidth[i]/2.0;
    }
    for(i=nbDimensions-4;i<nbDimensions;i++) {
      temp=point1->valeursTraitees[i]-point2->valeursTraitees[i];
      temp*=temp;
      valeur1+=temp*point1->inverseBandwidth[i];
      valeur2+=temp*point2->inverseBandwidth[i];
    }
    valeur=0.5*exp(-valeur1)+0.5*exp(-valeur2);
    if(valeur>valeurMax) {
      valeurMax=valeur;
    }
  }
  valeur=valeurMax;
  return valeur;
}*/
/*noyau special MNIST (bandwidth de 37 pour distribution exponentielle)*/
/*double radialGaussianKernel(int nbDimensions,individu_t *point1,individu_t *point2) {
  int i,j;
  double temp,valeur;

  valeur=0;
  
  for(i=0;i<nbDimensions;i++) {
    temp=point1->valeursTraitees[i]-point2->valeursTraitees[i];
    temp*=temp;
    valeur+=temp*point1->inverseBandwidth[i];
  }

  int distance,rayon,nbValeurs;
  double valeurMoyenne1,valeurMoyenne2;
  for(rayon=2;rayon<11;rayon++) {
    nbValeurs=0;
    valeurMoyenne1=0;
    valeurMoyenne2=0;
    for(i=0;i<28;i++) {
      for(j=0;j<28;j++) {
	distance=(int)round(sqrt((double)((i-14)*(i-14)+(j-14)*(j-14))));
	if(distance==rayon) {
	  valeurMoyenne1+=point1->valeursTraitees[i*28+j];
	  valeurMoyenne2+=point2->valeursTraitees[i*28+j];
	  nbValeurs++;
	}
      }
    }
    temp=(valeurMoyenne1-valeurMoyenne2)/(double)nbValeurs;
    valeur+=20*temp*temp*point1->inverseBandwidth[0];
  }
 
  int ligne;
  for(ligne=0;ligne<28;ligne++) {
    valeurMoyenne1=0;
    valeurMoyenne2=0;
    for(i=0;i<28;i++) {
      valeurMoyenne1+=point1->valeursTraitees[ligne*28+j];
      valeurMoyenne2+=point2->valeursTraitees[ligne*28+j];
    }
    temp=(valeurMoyenne1-valeurMoyenne2)/28.0;
    valeur+=10*temp*temp*point1->inverseBandwidth[0];
  }
  
  int colonne;
  for(colonne=0;colonne<28;colonne++) {
    valeurMoyenne1=0;
    valeurMoyenne2=0;
    for(i=0;i<28;i++) {
      valeurMoyenne1+=point1->valeursTraitees[i*28+colonne];
      valeurMoyenne2+=point2->valeursTraitees[i*28+colonne];
    }
    temp=(valeurMoyenne1-valeurMoyenne2)/28.0;
    valeur+=10*temp*temp*point1->inverseBandwidth[0];
    }
  valeur=exp(-valeur);
  
  return valeur;
  }*/
/*noyau normal*/
double radialGaussianKernel(int nbDimensions,individu_t *point1,individu_t *point2) {
  double valeur;
  int i;
  double temp,valeur1,valeur2;
  valeur1=0;
  valeur2=0;
  
  for(i=0;i<nbDimensions;i++) {
    temp=point1->valeursTraitees[i]-point2->valeursTraitees[i];
    temp*=temp;
    valeur1+=temp*point1->inverseBandwidth[i];
    valeur2+=temp*point2->inverseBandwidth[i];
  }
  valeur=0.5*exp(-valeur1)+0.5*exp(-valeur2);
  return valeur;
}

/******************************************************************/
/*                                                                */
/*Procedure de calcul du noyau d'Epanenchnikov                    */
/*bandwidth contient les variances pour les differentes dimensions*/
/*                                                                */
/******************************************************************/

double EpanechnikovKernel(int nbDimensions,individu_t *point1,individu_t *point2) {
  double valeur1,valeur2;
  int i;
  double temp,distanceCarre;

  distanceCarre=0;
  for(i=0;i<nbDimensions;i++) {
    temp=point1->valeursTraitees[i]-point2->valeursTraitees[i];
    distanceCarre+=temp*temp;
  }
  if(distanceCarre<1.0/point1->inverseBandwidth[0]) {
    valeur1=1.0-distanceCarre*point1->inverseBandwidth[0];
  } else {
    valeur1=0;
  }
  if(distanceCarre<1.0/point1->inverseBandwidth[0]) {
    valeur2=1.0-distanceCarre*point1->inverseBandwidth[0];
  } else {
    valeur2=0;
  }
  
  return 0.5*valeur1+0.5*valeur2;
}

/******************************************************************/
/*                                                                */
/*Procedure de calcul du noyau polynomial                         */
/*bandwidth contient les variances pour les differentes dimensions*/
/*                                                                */
/******************************************************************/

double polynomialKernel(int nbDimensions,individu_t *point1,individu_t *point2,int choixNoyau) {
  double valeur;
  int i;
  double temp,valeur1,valeur2;
  valeur1=0;
  valeur2=0;
  for(i=0;i<nbDimensions;i++) {
    temp=point1->valeursTraitees[i]*point2->valeursTraitees[i];
    valeur1+=temp*point1->inverseBandwidth[i];
    valeur2+=temp*point2->inverseBandwidth[i];
  }
  switch (choixNoyau) {
  case POLYNOMIAL2 :
    valeur=0.5*pow(valeur1,2.0)+0.5*pow(valeur2,2.0);
    break;
  case POLYNOMIAL3 :
    valeur=0.5*pow(valeur1,3.0)+0.5*pow(valeur2,3.0);
    break;
  case POLYNOMIAL4 : 
    valeur=0.5*pow(valeur1,4.0)+0.5*pow(valeur2,4.0);
    break;
  }
  return valeur;
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul du bandwidth par estimation de la densite*/
/*                                                             */
/***************************************************************/

void calculFixedBandwidth(int nbDimensions,int nbIndividus,individu_t **individus,double bandwidth) {
  int i,j;
  double temp;
  temp=1.0/(bandwidth*bandwidth);
  for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbDimensions;j++) {
      individus[i]->inverseBandwidth[j]=temp;
    }
  }
}

/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculAdaptiveBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet,
			     int kPlusProches,double scale)  {
  int i,j,k;
  double temp,distance,kiemeDistance,*logfi,logG;
  individu_t **individusCopie;

  /*allocation memoire*/
  logfi=(double *)ckalloc(sizeof(double)*tailleTrainingSet);
  individusCopie=(individu_t **)ckalloc(sizeof(individu_t *)*tailleTrainingSet);
  /*fin allocation memoire*/
  
  for(i=0;i<tailleTrainingSet;i++) {
    individusCopie[i]=trainingSet[i];
  }
  logG=0;
  for(i=0;i<tailleTrainingSet;i++) {
    for(j=0;j<tailleTrainingSet;j++) {
      distance=0;
      for(k=0;k<nbDimensions;k++) {
	temp=trainingSet[i]->valeursTraitees[k]-individusCopie[j]->valeursTraitees[k];
	distance+=temp*temp;
      }
      individusCopie[j]->valeurTemp=sqrt(distance);
    }
    qsort(individusCopie,tailleTrainingSet,sizeof(individu_t *),compareIndividus);
    kiemeDistance=individusCopie[kPlusProches-1]->valeurTemp;
    logfi[i]=(double)nbDimensions*log(kiemeDistance);
    logG+=logfi[i];
  }
  logG/=(double)tailleTrainingSet;
  
  for(i=0;i<tailleTrainingSet;i++) {
    temp=exp(logfi[i]-logG)*scale;
    for(j=0;j<nbDimensions;j++) {
      trainingSet[i]->inverseBandwidth[j]=1.0/(temp*temp);
    }
  }

  /*des allocation memoire*/
  ckfree(logfi);
  ckfree(individusCopie);
  /*fin desallocation memoire*/
}

/*************************************************/
/*                                               */
/*Procedure de calcul du bandwidth Component-Wise*/
/*                                               */
/*************************************************/

void calculAdaptiveBandwidthComponentWise(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet,
				   	  int kPlusProches,double scale)  {
  int i,j,dimension;
  double temp,distance,**logfij,*logG;
  individu_t **individusCopie;

  /*allocation memoire*/
  logfij=(double **)ckalloc(sizeof(double *)*tailleTrainingSet);
  for(i=0;i<tailleTrainingSet;i++) {
  	logfij[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
  }
  logG=(double *)ckalloc(sizeof(double)*nbDimensions);
  individusCopie=(individu_t **)ckalloc(sizeof(individu_t *)*tailleTrainingSet);
  /*fin allocation memoire*/
  

  for(dimension=0;dimension<nbDimensions;dimension++) {
  	logG[dimension]=0;
  	for(i=0;i<tailleTrainingSet;i++) {
  	        individusCopie[i]=trainingSet[i];
  	}
  	for(i=0;i<tailleTrainingSet;i++) {
 	   for(j=0;j<tailleTrainingSet;j++) {
		individusCopie[j]->valeurTemp=fabs(trainingSet[i]->valeursTraitees[dimension]-individusCopie[j]->valeursTraitees[dimension]);
    	    }
 	    qsort(individusCopie,tailleTrainingSet,sizeof(individu_t *),compareIndividus);
            logfij[i][dimension]=(double)nbDimensions*individusCopie[kPlusProches-1]->valeurTemp;
    	    logG[dimension]+=logfij[i][dimension];
        }
  }
  for(j=0;j<nbDimensions;j++) {
	logG[j]/=(double)tailleTrainingSet;
  }  
  for(i=0;i<tailleTrainingSet;i++) {
    for(j=0;j<nbDimensions;j++) {
      temp=exp(logfij[i][j]-logG[j])*scale;
      trainingSet[i]->inverseBandwidth[j]=1.0/(temp*temp);
    }
  }

  /*des allocation memoire*/
  for(i=0;i<tailleTrainingSet;i++) {
  	ckfree(logfij[i]); 
  }
  ckfree(logfij);
  ckfree(logG);
  ckfree(individusCopie);
  /*fin desallocation memoire*/
}

/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculMMDensityBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet)  {
  int i,j;
  double temp,*logfi,logG;
  
  /*allocation memoire*/
  logfi=(double *)ckalloc(sizeof(double)*tailleTrainingSet);
  /*fin allocation memoire*/

  logG=0;
  for(i=0;i<tailleTrainingSet;i++) {
      logfi[i]=log(trainingSet[i]->vraisemblance);
      logG+=logfi[i];
  }
  logG/=(double)tailleTrainingSet;
  
  for(i=0;i<tailleTrainingSet;i++) {
    temp=exp(-0.5*(logfi[i]-logG));
    for(j=0;j<nbDimensions;j++) {
      trainingSet[i]->inverseBandwidth[j]=1.0/(temp*temp);
    }
  }

  /*des allocation memoire*/
  ckfree(logfi);
  /*fin desallocation memoire*/

}

/**********************************************************************/
/*                                                                    */
/*Bandwidth adaptatif d'Abramson avec une densite pilote a noyau      */
/*(utilise par -cm kkmeans)                                           */
/*                                                                    */
/**********************************************************************/
/*2026-10-08 : calculMMDensityBandwidth donnait h_i = (f_i/G)^-1/2 sans
  bandwidth global (h=1 dans l'unite des donnees, quelle que soit leur
  echelle), et la densite f_i, estimee par modele de melange, ignorait la
  composante bruit : les points d'un nuage absorbe par le bruit avaient
  f_i ~ 1e-100, d'ou des inverseBandwidth de 1e-100 a 1e+55 et un noyau
  nul entre points distincts (kkmeans : groupes vides, inertie NaN).
  Ici : h0^2 = mediane des distances au carre entre individus / 2 (le noyau
  vaut exp(-1) pour une paire a distance mediane), densite pilote
  f_i = moyenne_j exp(-d_ij^2/(2 h0^2)) (toujours >= 1/n), puis
  h_i = h0 (f_i/G)^-1/2, G moyenne geometrique des f_i.
  radialGaussianKernel calcule exp(-sum d^2 * inverseBandwidth) :
  inverseBandwidth = 1/(2 h_i^2).*/

static double distanceCarre(int nbDimensions,individu_t *point1,individu_t *point2) {
  int k;
  double temp,somme=0;
  for(k=0;k<nbDimensions;k++) {
    temp=point1->valeursTraitees[k]-point2->valeursTraitees[k];
    somme+=temp*temp;
  }
  return somme;
}

static int compareDoubles(const void *a,const void *b) {
  double x=*(const double *)a,y=*(const double *)b;
  return (x<y)?-1:((x>y)?1:0);
}

void calculPilotAdaptiveBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet) {
  int i,j,k,nbPaires;
  double temp,h02,logG,*distances2,*logfi;

  nbPaires=tailleTrainingSet*(tailleTrainingSet-1)/2;

  /*allocation memoire*/
  distances2=(double *)ckalloc(sizeof(double)*(nbPaires>0?nbPaires:1));
  logfi=(double *)ckalloc(sizeof(double)*tailleTrainingSet);
  /*fin allocation memoire*/

  /*bandwidth global : mediane des distances au carre*/
  k=0;
  for(i=0;i<tailleTrainingSet;i++) {
    for(j=i+1;j<tailleTrainingSet;j++) {
      distances2[k++]=distanceCarre(nbDimensions,trainingSet[i],trainingSet[j]);
    }
  }
  h02=0;
  if(nbPaires>0) {
    qsort(distances2,nbPaires,sizeof(double),compareDoubles);
    h02=distances2[nbPaires/2]/2.0;
    if(h02<=0) {
      /*plus de la moitie des individus confondus : distance moyenne*/
      for(k=0;k<nbPaires;k++) {
	h02+=distances2[k];
      }
      h02/=2.0*(double)nbPaires;
    }
  }
  if(h02<=0) {
    /*tous les individus sont identiques*/
    h02=1.0;
  }

  /*densite pilote (noyau gaussien de bandwidth h0) et moyenne geometrique*/
  logG=0;
  for(i=0;i<tailleTrainingSet;i++) {
    temp=0;
    for(j=0;j<tailleTrainingSet;j++) {
      temp+=exp(-distanceCarre(nbDimensions,trainingSet[i],trainingSet[j])/(2.0*h02));
    }
    logfi[i]=log(temp/(double)tailleTrainingSet);
    logG+=logfi[i];
  }
  logG/=(double)tailleTrainingSet;

  /*h_i^2 = h0^2 * G/f_i*/
  for(i=0;i<tailleTrainingSet;i++) {
    temp=h02*exp(logG-logfi[i]);
    for(j=0;j<nbDimensions;j++) {
      trainingSet[i]->inverseBandwidth[j]=1.0/(2.0*temp);
    }
  }

  /*desallocation memoire*/
  ckfree(distances2);
  ckfree(logfi);
  /*fin desallocation memoire*/
}

/**********************************************************************/
/*                                                                    */
/*Procedure de calcul du bandwidth (ordre de grandeur de l'ecart-type)*/
/*                                                                    */
/**********************************************************************/

void calculNNAdaptiveBandwidth(int nbDimensions,int tailleTrainingSet,individu_t **trainingSet,
			       int kPlusProches,double scale,int accelerationCalculs,int kPlusProchesVoisins)  {
  int i,j,k;
  double temp,distance,kiemeDistance,*min,*max;
  individu_t **individusCopie;
  
  if((kPlusProches>kPlusProchesVoisins)&&(accelerationCalculs==OUI)) {
    printf("boost too large !\n");
    exit(0);
  }
  
  /*allocation memoire*/
  min=(double *)ckalloc(sizeof(double)*nbDimensions);
  max=(double *)ckalloc(sizeof(double)*nbDimensions);
  individusCopie=(individu_t **)ckalloc(sizeof(individu_t *)*tailleTrainingSet);
  /*fin allocation memoire*/

  for(i=0;i<tailleTrainingSet;i++) {
    individusCopie[i]=trainingSet[i];
  }

  for(i=0;i<tailleTrainingSet;i++) { 
    if(accelerationCalculs==NON) {
      for(j=0;j<tailleTrainingSet;j++) {
	distance=0;
	for(k=0;k<nbDimensions;k++) {
	  temp=trainingSet[i]->valeursTraitees[k]-individusCopie[j]->valeursTraitees[k];
	  distance+=temp*temp;
	}
	distance=sqrt(distance);
	individusCopie[j]->valeurTemp=distance;
      }      
      qsort(individusCopie,tailleTrainingSet,sizeof(individu_t *),compareIndividus);
      kiemeDistance=individusCopie[kPlusProches]->valeurTemp;
    } else {
      kiemeDistance=0;
      for(k=0;k<nbDimensions;k++) {
	temp=trainingSet[i]->kNearestNeighbours[kPlusProches]->valeursTraitees[k]-
	  trainingSet[i]->valeursTraitees[k];
	kiemeDistance=temp*temp;
      }    
    }

    for(j=0;j<nbDimensions;j++) {
      min[j]=DBL_MAX;
      max[j]=-DBL_MAX;
    }
    for(j=0;j<nbDimensions;j++) {
      for(k=0;k<=kPlusProches;k++) {
	if(individusCopie[k]->valeursTraitees[j]>max[j]) {
	  max[j]=individusCopie[k]->valeursTraitees[j];
	} else if(individusCopie[k]->valeursTraitees[j]<min[j]) {
	  min[j]=individusCopie[k]->valeursTraitees[j];
	}
      }
      if(max[j]-min[j]>ZERO_LIMIT) {
	trainingSet[i]->inverseBandwidth[j]=1.0/pow(max[j]-min[j],2.0)/(scale*scale);
      } else {
	trainingSet[i]->inverseBandwidth[j]=1.0/pow(sqrt((double)nbDimensions)/kiemeDistance,2.0)/(double)nbDimensions/(scale*scale);
      }
      //trainingSet[i]->inverseBandwidth[j]=1.0/pow(max[j],2.0)/(double)nbDimensions;
      //trainingSet[i]->inverseBandwidth[j]=pow(sqrt((double)nbDimensions)/kiemeDistance,2.0);
    }
  } 

  /*desallocation memoire*/
  ckfree(min);
  ckfree(max);
  ckfree(individusCopie);
  /*fin desallocation memoire*/
}
