#include "mainClusPack.h"
#include "space.h"
#include "tcl.h"

/********************************************************************************************/
/*                                                                                          */
/*Procedure de recherche des k plus proches voisins pour un cluster avec stockage des pivots*/
/*                                                                                          */
/********************************************************************************************/

void kNearestNeighbours(int nbPivots,
			int kPlusProches,
			int nbIndividus,
			int nbDimensions,
			individu_t **individus,
			cluster_t *cluster) 
{
  int i,j,k;
  individu_t **individusCopie;
  double *distancesMax;

  printf("kplusproches : %d\n",kPlusProches);
  //  getchar();
  /*allocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    individus[i]->kNearestNeighbours=(individu_t **)ckalloc(sizeof(individu_t *)*kPlusProches);
    individus[i]->distancesPivots=(double *)ckalloc(sizeof(double)*nbPivots);
  }
  cluster->pivots=(individu_t **)ckalloc(sizeof(individu_t *)*nbPivots);
  distancesMax=(double *)ckalloc(sizeof(double)*nbIndividus);
  individusCopie=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  /*fin allocation memoire*/

  /*selection aleatoire de nbPivots*/
  for(i=0;i<nbPivots;i++) {
    cluster->pivots[i]=individus[(int)floor((double)nbIndividus*runif(0,1))];
  }
  double temp,base;
  
  /*calcul des noyaux entre les pivots et les individus*/
  //  min=DBL_MAX;
  for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbPivots;j++) {
      //individusCopie[j]=cluster->pivots[j];
      individus[i]->distancesPivots[j]=0;
      for(k=0;k<nbDimensions;k++) { 
	temp=individus[i]->valeursTraitees[k]-cluster->pivots[j]->valeursTraitees[k];
	individus[i]->distancesPivots[j]+=temp*temp;
      }
      individus[i]->distancesPivots[j]=sqrt(individus[i]->distancesPivots[j]);
      //printf("dist : %d %f\n",j,distancesIndividusPivots[i][j]);
      /*if((0<distancesIndividusPivots[i][j])&&(distancesIndividusPivots[i][j]<min)) {
	min=distancesIndividusPivots[i][j];
	}
	individusCopie[j]->valeurTemp=distancesIndividusPivots[i][j];*/
    }
    //qsort(individusCopie,nbPivots,sizeof(individu_t *),compareIndividus);
    //distancesMax[i]=individusCopie[kPlusProches-1]->valeurTemp;
    //distancesMax[i]=min;
    //printf("min : %f\n",min);
    //getchar();
  }
  //printf("min : %f\n",min);
  base=1.0;
  //getchar();
  /*calcul d'une borne sup pour le kieme voisin de chaque individus*/
  int compteur,individuSelectionne;
  double facteur;
  for(i=0;i<nbIndividus;i++) {
    facteur=1.0;
  debutBouclekNearestNeighbours1:
    distancesMax[i]=base*facteur;
    compteur=0;
    for(j=0;j<nbIndividus;j++) {
      individuSelectionne=OUI;
      for(k=0;k<nbPivots;k++) {
	//printf("%f\n",fabs(distancesIndividusPivots[i][k]-distancesIndividusPivots[j][k]));
	
	if(fabs(individus[i]->distancesPivots[k]-individus[j]->distancesPivots[k])>distancesMax[i]) {
	  individuSelectionne=NON;
	  break;
	}
      }
      if(individuSelectionne==OUI) {
	individusCopie[compteur]=individus[j];
	compteur++;
      }
      //printf("compteur : %d\n",compteur);
      //getchar();
    }//printf("compteur avant : %d\n",compteur);
    if(compteur<kPlusProches) {
      facteur*=2.0;
      goto debutBouclekNearestNeighbours1;
    }
    //printf("compteur : %d\n",compteur);
      //getchar();
    for(j=0;j<compteur;j++) {
      individusCopie[j]->valeurTemp=0;
      for(k=0;k<nbDimensions;k++) {
	temp=individus[i]->valeursTraitees[k]-individusCopie[j]->valeursTraitees[k];
	individusCopie[j]->valeurTemp+=temp*temp;
      } 
    }
    qsort(individusCopie,compteur,sizeof(individu_t *),compareIndividus);
    for(j=0;j<kPlusProches;j++) {
      individus[i]->kNearestNeighbours[j]=individusCopie[j];
    }
  }

  /*desallocation memoire*/
  ckfree(distancesMax);
  ckfree(individusCopie);
  /*fin desallocation memoire*/
  // getchar();
}

/*****************************************************************************************/
/*                                                                                       */
/*Procedure de recherche des k plus proches voisins pour un individu dans un cluster*/
/*                                                                                       */
/*****************************************************************************************/

void kNearestNeighboursInOtherCluster(int nbPivots,int kPlusProches,int nbDimensions,individu_t *individu,
				      cluster_t *cluster) {
  int j,k,nbIndividus;
  individu_t **individusCopie;
  double *distancesIndividuPivots,distanceMax,base;

  nbIndividus=cluster->nbIndividus;

  //printf("kplusproches dans cluster : %d\n",kPlusProches);

  /*allocation memoire*/
  individu->kNearestNeighbours=(individu_t **)ckalloc(sizeof(individu_t *)*kPlusProches);
  distancesIndividuPivots=(double *)ckalloc(sizeof(double)*nbPivots);
  individusCopie=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  /*fin allocation memoire*/

  double temp;
  /*calcul des noyaux entre les pivots et les individus*/
  for(j=0;j<nbPivots;j++) {
    //individusCopie[j]=cluster->pivots[j];
    distancesIndividuPivots[j]=0;
    for(k=0;k<nbDimensions;k++) { 
      temp=individu->valeursTraitees[k]-cluster->pivots[j]->valeursTraitees[k];
      distancesIndividuPivots[j]+=temp*temp;
    } 
    distancesIndividuPivots[j]=sqrt(distancesIndividuPivots[j]);
    //individusCopie[j]->valeurTemp=distancesIndividuPivots[j];
  }
  //qsort(individusCopie,nbPivots,sizeof(individu_t *),compareIndividus);
  //  distanceMax=individusCopie[kPlusProches-1]->valeurTemp;
  base=1.0;

  /*calcul d'une borne sup pour le kieme voisin de chaque individus*/
  int compteur,individuSelectionne;
  double facteur=1.0;
  
 debutBouclekNearestNeighbours2:
  distanceMax=base*facteur;
  compteur=0;
  for(j=0;j<nbIndividus;j++) {
    individuSelectionne=OUI;
    for(k=0;k<nbPivots;k++) {
      if(fabs(distancesIndividuPivots[k]-cluster->individus[j]->distancesPivots[k])>distanceMax) {
	individuSelectionne=NON;
	break;
      }
    }
    if(individuSelectionne==OUI) {
      individusCopie[compteur]=cluster->individus[j];
      compteur++;
    }
  }
  
  if(compteur<kPlusProches) {
    facteur*=2.0;
    goto debutBouclekNearestNeighbours2;
  }
  for(j=0;j<compteur;j++) {
    individusCopie[j]->valeurTemp=0;
    for(k=0;k<nbDimensions;k++) {
      temp=individu->valeursTraitees[k]-individusCopie[j]->valeursTraitees[k];
      individusCopie[j]->valeurTemp+=temp*temp;
    } 
  }
  //printf("compteur : %d\n",compteur);
  qsort(individusCopie,compteur,sizeof(individu_t *),compareIndividus);
  for(j=0;j<kPlusProches;j++) {
    individu->kNearestNeighbours[j]=individusCopie[j];
  }

  /*desallocation memoire*/
  ckfree(distancesIndividuPivots);
  ckfree(individusCopie);
  /*fin desallocation memoire*/
}
