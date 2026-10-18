#include "mainClusPack.h"
#include "mixtureModels.h"
#include "matrice.h"
#include "symeig.h"
#include "emGaussian.h"
#include "eigenMixtureModels.h"

/*************************************************/
/*                                               */
/*Clustering par eigen mixture models            */
/*Etape 1 : on surestime le nombre de composantes*/
/*Etape 2 : on fusionne les composantes par PCA  */
/*                                               */
/*************************************************/

void eigenMixtureModels(int nbIndividus,int nbDimensions,individu_t **individus,int nbIterationsMax,
			int nbClustersMax,int *nbClusters,int diag) {
  int i,j,typeDonnees,nbEssais,nbEssaisMax;
  mixtureModel_t mixtureModel;

  nbEssaisMax=20;
  /*if(nbClustersMax==0) {
    nbClustersMax=(int)sqrt((double)nbIndividus);
    if((nbClustersMax<nbIndividus/(nbDimensions+1))&&(nbIndividus/(nbDimensions+1)<60))
      {
	nbClustersMax=nbIndividus/(nbDimensions+1);
      }
  } else {
    printf("nbClustersMax : %d\n",nbClustersMax);
    }*/
  
  // mixtureModel.nbClusters=nbClustersMax;

  typeDonnees=GAUSSIAN;

  //mixtureModel.nbClusters=nbClustersMax;

  /*on effectue un clustering par nbClustersMax modeles de melanges*/
  nbEssais=0;
  mixtureModel.destroy=0;
  while(nbEssais<nbEssaisMax) {
    if(clusteringMixtureModels(nbIndividus,nbDimensions,individus,&mixtureModel,nbIterationsMax,AIC,
			       typeDonnees,diag,-1,OUI)==SUCCES) {
      break;
    }
    nbEssais++;
  }

  /*creation des clusters*/

  printf("fin du modele de melanges\n");
  fflush(stdout);

  /*calcul des vraisemblances pour chaque individu*/
  calculVraisemblancesPointsGaussienne(nbIndividus,nbDimensions,individus,
				       mixtureModel.nbClusters-1,mixtureModel.clusters,diag);

  /*calcul des nouvelles coordonnees*/
  int nbDimensions2;
  double **points,**matriceCovariance;
  
  nbDimensions2=mixtureModel.nbClusters-1;
  printf("nbDimensions2 : %d\n",nbDimensions2);
  /*allocation memoire*/
  points=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    points[i]=(double *)ckalloc(sizeof(double)*nbDimensions2);
  }
  matriceCovariance=(double **)ckalloc(sizeof(double *)*nbDimensions2);
  for(i=0;i<nbDimensions2;i++) {
    matriceCovariance[i]=(double *)ckalloc(sizeof(double)*nbDimensions2);
  }
  /*fin allocation memoire*/

  double somme;
  for(i=0;i<nbIndividus;i++) {
    somme=0;
    for(j=0;j<nbDimensions2;j++) {
      somme+=exp(individus[i]->vrais[j]);
    }
    for(j=0;j<nbDimensions2;j++) {
      points[i][j]=exp(individus[i]->vrais[j]);
    }
  }

  /*calcul de la matrice de covariance*/
  calculMatriceCovariance(nbDimensions2,nbIndividus,points,matriceCovariance);
  double *valeursPropres,**vecteursPropres;
  
  /*allocation memoire*/
  valeursPropres=(double *)ckalloc(sizeof(double)*nbDimensions2);
  vecteursPropres=(double **)ckalloc(sizeof(double *)*nbDimensions2);
  for(i=0;i<nbDimensions2;i++) {
    vecteursPropres[i]=(double *)ckalloc(sizeof(double)*nbDimensions2);
  }
  /*fin allocation memoire*/

  /*calcul des vecteurs propres*/
  symeig(matriceCovariance,nbDimensions2,vecteursPropres,valeursPropres);
  printf("matrice de covariance\n");
  for(i=0;i<nbDimensions2;i++) {
    for(j=0;j<nbDimensions2;j++) {
      printf("%f ",matriceCovariance[i][j]);
    }
    printf("\n");
  }
  for(i=0;i<nbDimensions2;i++) {
    printf("eigenvalues : %f\n",valeursPropres[i]);
  }


  /*printf("matrice de covariance\n");
  for(i=0;i<nbDimensions2;i++) {
    for(j=0;j<nbDimensions2;j++) {
      printf("%.2f ",matriceCovariance[i][j]);
    }
    printf("\n");
  }

  for(i=0;i<nbDimensions2;i++) {
    printf("eigenvalues : %lf\n",valeursPropres[i]);
    }*/
  for(i=0;i<nbDimensions2;i++) {
    for(j=0;j<nbDimensions2;j++) {
      printf("%f ",vecteursPropres[i][j]);
    }
    printf("\n");
  }  

  int *vecteursDesComposantes;

  /*allocation memoire*/
  vecteursDesComposantes=(int *)ckalloc(sizeof(int)*nbDimensions2);
  /*fin allocation memoire*/

  /*affectation des composantes aux differents vecteurs propres*/
  affectationComposantesMelangeAuxVecteursPropres(nbDimensions2,vecteursPropres,vecteursDesComposantes);

  /*for(i=0;i<nbDimensions2;i++) {
    printf("dimension %d affectee au vecteur : %d\n",i,vecteursDesComposantes[i]);
    }*/
  
  /*creation des clusters*/
  for(i=0;i<mixtureModel.nbClusters-1;i++) {
    for(j=0;j<mixtureModel.clusters[i].nbIndividus;j++) {
      if(vecteursDesComposantes[i]==NOT_DEFINED) {
	mixtureModel.clusters[i].individus[j]->cluster=i+1+mixtureModel.nbClusters;
      } else {
	mixtureModel.clusters[i].individus[j]->cluster=vecteursDesComposantes[i];
      }
    }
  }

  /*for(i=0;i<nbIndividus;i++) {
    printf("cluster de %d : %d\n",i,individus[i]->cluster);
    }*/


  /*a affiner plus tard*/
  *nbClusters=2*mixtureModel.nbClusters+1;

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    ckfree(points[i]);
  }
  ckfree(points);
  for(i=0;i<nbDimensions2;i++) {
    ckfree(vecteursPropres[i]);
  }
  ckfree(vecteursPropres);
  ckfree(valeursPropres);
  ckfree(vecteursDesComposantes);
  /*fin desallocation memoire*/
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul de la pvaleur d'avoir la somme des carres*/
/*de k coordonnees d'un vecteur propre normalise superieure a x*/
/*                                                             */
/***************************************************************/

double pvaleursKcoordonneesVecteurPropre(int nbDimensions,int k,double x) {
  int i,nbIterations,compteur;
  double Xk,Xp_k,t;

  t=sqrt(x)-sqrt((double)k/(double)nbDimensions)*0.5;
  printf("t : %f\n",t);
  /*if(t<0) {
    return facteur1=1.0;
  } else {
    facteur1=2.0*exp(-t*t*(double)nbDimensions/2.0);
    }*/
   
  nbIterations=100000;
  compteur=0;
  for(i=0;i<nbIterations;i++) {
    Xk=rchisq((double)k);
    Xp_k=rchisq((double)(nbDimensions-k));
    if(Xk>=x/(1.0-x)*Xp_k) {
      compteur++;
    }
  }
  double facteur2;
  facteur2=(double)nbDimensions*log((double)nbDimensions)-(double)k*log((double)k)-
    (double)(nbDimensions-k)*log((double)(nbDimensions-k))-(double)k*log(2.0);
  facteur2=exp(facteur2)*sqrt((double)nbDimensions/(2.0*M_PI*(double)(k*(nbDimensions-k))));
  //  printf("DIM : %d ; fact : %f\n",nbDimensions,facteur);
  //return (double)compteur/(double)nbIterations;
  /*if((double)k>=10.0*log((double)nbDimensions)) {
    return facteur1*facteur2;
    } else {*/
    return (double)compteur/(double)nbIterations*facteur2;
    //}
}

/*******************************************************/
/*                                                     */
/*Procedure de calcul de comparaison pour le classement*/
/*des valeurs propres de maniere croissante            */
/*                                                     */
/*******************************************************/

int compareCoordonnees(void const *a,void const *b) {
  individu_t **pa=(individu_t **)a;
  individu_t **pb=(individu_t **)b;
  double diff=(*pa)->valeurTemp-(*pb)->valeurTemp;

  if(diff<0) {
    return 1;
  } else {
    return 0;
  }
}

/****************************************************************/
/*                                                              */
/*Procedure d'affectation des differentes composantes du melange*/
/*aux differents vecteurs propres                               */
/*                                                              */
/****************************************************************/

void affectationComposantesMelangeAuxVecteursPropres(int nbDimensions,double **vecteursPropres,
						     int *vecteursDesComposantes) {
  int i,uneComposanteSansVecteur;
  individu_t **coordonneesAClasser;
  uneComposanteSansVecteur=OUI;

  for(i=0;i<nbDimensions;i++) {
    vecteursDesComposantes[i]=NOT_DEFINED;
  }
  /*allocation memoire*/
  coordonneesAClasser=(individu_t **)ckalloc(sizeof(individu_t *)*nbDimensions);
  for(i=0;i<nbDimensions;i++) {
    coordonneesAClasser[i]=(individu_t *)ckalloc(sizeof(individu_t)*1);
  }
  /*fin allocation memoire*/
  
  int meilleurNbValeursAjoutees,nbValeursAjoutees,nbDimensionsActualise,vecteurCourant;
  double sommeCarres,pvaleurMin,pvaleur;

  /*le vecteur avec la plus grande valeur propre est le dernier, 
    c'est donc par celui-la que l'on commence*/
  vecteurCourant=nbDimensions-1;
  while(vecteurCourant>=0) {
    uneComposanteSansVecteur=NON;
    for(i=0;i<nbDimensions;i++) {
      if(vecteursDesComposantes[i]==NOT_DEFINED) {
	uneComposanteSansVecteur=OUI;
	break;
      }
    }
    if(uneComposanteSansVecteur==NON) {
      break;
    }

    /*on normalise les coordonnees restantes*/
    sommeCarres=0;
    for(i=0;i<nbDimensions;i++) {
      if(vecteursDesComposantes[i]==NOT_DEFINED) {
	sommeCarres+=vecteursPropres[i][vecteurCourant]*vecteursPropres[i][vecteurCourant];
      }
    }
    sommeCarres=sqrt(sommeCarres);
    for(i=0;i<nbDimensions;i++) {
      if(vecteursDesComposantes[i]==NOT_DEFINED) {
	vecteursPropres[i][vecteurCourant]/=sommeCarres;
      }
    }

    nbDimensionsActualise=0;
    for(i=0;i<nbDimensions;i++) {
      if(vecteursDesComposantes[i]==NOT_DEFINED) {
	coordonneesAClasser[nbDimensionsActualise]->id=i;
	/*modif*/
	coordonneesAClasser[nbDimensionsActualise]->valeurTemp=vecteursPropres[i][vecteurCourant];
	nbDimensionsActualise++;
      }
    }
    /*classement des coordonnees*/
    qsort(coordonneesAClasser,nbDimensionsActualise,sizeof(individu_t *),compareCoordonnees);

    sommeCarres=0;
    nbValeursAjoutees=0;
    pvaleurMin=1.0;
    meilleurNbValeursAjoutees=-1;
    for(i=0;i<nbDimensionsActualise;i++) {
      if(coordonneesAClasser[i]->valeurTemp<=0) {
	break;
      }
      nbValeursAjoutees++;
      sommeCarres+=coordonneesAClasser[i]->valeurTemp*coordonneesAClasser[i]->valeurTemp;
      
      pvaleur=pvaleursKcoordonneesVecteurPropre(nbDimensionsActualise,nbValeursAjoutees,sommeCarres);
      printf("%f (pvaleur : %f)\n",coordonneesAClasser[i]->valeurTemp,pvaleur);
      if(pvaleur<=pvaleurMin) {
	pvaleurMin=pvaleur;
	meilleurNbValeursAjoutees=nbValeursAjoutees;
      }
    }
    if(pvaleurMin>0.01) {
      break;
    }
    /*mise a jour*/
    meilleurNbValeursAjoutees=nbValeursAjoutees;
    printf("meilleur nb valeurs ajoutees : %d (pvaleur : %.4f)\n",meilleurNbValeursAjoutees,pvaleurMin);
    for(i=0;i<meilleurNbValeursAjoutees;i++) {
      vecteursDesComposantes[coordonneesAClasser[i]->id]=vecteurCourant;
    }
    vecteurCourant--;
  }
  
  /*allocation memoire*/
  for(i=0;i<nbDimensions;i++) {
    ckfree(coordonneesAClasser[i]);
  }
  ckfree(coordonneesAClasser);
  /*fin allocation memoire*/
}
