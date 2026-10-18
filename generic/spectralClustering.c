#include "mainClusPack.h"
#include "symeig.h"
#include "mixtureModels.h"
#include "matrice.h"
#include "kernels.h"
#include "densityEstimation.h"

/*************************************************/
/*                                               */
/*Procedure de realisation de spectral clustering*/
/*en partant d'un debut de clustering            */ 
/*                                               */
/*************************************************/

void spectralClustering(int nbIndividus,int nbDimensions,individu_t **individus,int nbClustersSelection,
			int *nbClusters,int bandwidthType,double bandwidth,int diag) {
  int i,j,k,nbDimensionsSpectrales,kPlusProches,nbDimensionsSpectralesMax;
  double temp,**laplacien,scale;
  
  nbDimensionsSpectralesMax=10;
  kPlusProches=10;
  scale=1.0;
  if(nbDimensions<nbDimensionsSpectralesMax) {
    nbDimensionsSpectrales=nbDimensions;
  } else {
    nbDimensionsSpectrales=nbDimensionsSpectralesMax;
  }
  
  /*calcul du bandwidth*/
  if(bandwidthType==FIXED_BANDWIDTH) {
    calculFixedBandwidth(nbDimensions,nbIndividus,individus,bandwidth);
  } else if(bandwidthType==AB) {
    calculAdaptiveBandwidth(nbDimensions,nbIndividus,individus,kPlusProches,scale);
  } else if(bandwidthType==ABMM) {
    mixtureModel_t mixtureModel1;
    mixtureModelDensityEstimation(nbIndividus,nbDimensions,individus,diag,&mixtureModel1) ;
    printf("fin mm\n");
    calculMMDensityBandwidth(nbDimensions,nbIndividus,individus);
    printf("fin calcul densite\n");
  } 

  /*allocation memoire*/
  laplacien=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    laplacien[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
  }
  /*fin allocation memoire*/

  /*calcul de la matrice de distance et de la plus grande distance*/
  for(i=0;i<nbIndividus;i++) {
    for(j=i+1;j<nbIndividus;j++) {
      laplacien[i][j]=0;
      for(k=0;k<nbDimensions;k++) {
	temp=individus[i]->valeursTraitees[k]-individus[j]->valeursTraitees[k];
	laplacien[i][j]+=temp*temp;
      }
      laplacien[i][j]=kernel(nbDimensions,individus[i],individus[j],RADIAL_GAUSSIAN);
      laplacien[j][i]=laplacien[i][j];
    }
  }

  double somme;
  for(i=0;i<nbIndividus;i++) { 
    somme=0;
    for(j=0;j<nbIndividus;j++) {
      if(j!=i) {
	somme+=laplacien[i][j];
      }
    }
    laplacien[i][i]=-somme;
  }

  double *valeursPropres,**vecteursPropres;

  /*allocation memoire*/
  valeursPropres=(double *)ckalloc(sizeof(double)*nbIndividus);
  vecteursPropres=(double **)ckalloc(sizeof(double *)*(nbDimensionsSpectrales+1));
  for(i=0;i<nbDimensionsSpectrales+1;i++) {
    vecteursPropres[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
  }
  /*fin allocation memoire*/

  int dimension,iteration,nbIterations=10;
  double **A,**B,*vecteur1,*vecteur2;
  
  /*allocation memoire*/
  A=(double **)ckalloc(sizeof(double *)*nbIndividus);
  B=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    A[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
    B[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
  }
  vecteur1=(double *)ckalloc(sizeof(double)*nbIndividus);
  vecteur2=(double *)ckalloc(sizeof(double)*nbIndividus);
  /*fin allocation memoire*/
   
  /*calcul des vecteurs propres, le premier etant le vecteur (1,1,1,...,1)*/
  /*pour cela on effectue un spectral shift puis inverse la matrice et calcul les plus grandes valeurs propres*/
   
  /*inversion de la matrice*/
  for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbIndividus;j++) {
      A[i][j]=laplacien[i][j];
    }
    /*spectral shift*/
    A[i][i]+=1.0;
  }
  
  inversionMatrice(nbIndividus,A);
  /*for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbIndividus;j++) {
      printf("%.4f ",A[i][j]);
    }
    printf("\n");
    }*/
  /*nbIndividus=4;
  A[0][0]=6.0;
  A[0][1]=4.0;
  A[0][2]=1.0;
  A[0][3]=1.0;
  A[1][0]=4.0;
  A[1][1]=6.0;
  A[1][2]=1.0;
  A[1][3]=1.0;
  A[2][0]=1.0;
  A[2][1]=1.0;
  A[2][2]=5.0;
  A[2][3]=2.0;
  A[3][0]=1.0;
  A[3][1]=1.0;
  A[3][2]=2.0;
  A[3][3]=5.0;*/
  double denominateur;

  for(dimension=0;dimension<nbDimensionsSpectrales+1;dimension++) {
    for(j=0;j<nbIndividus;j++) {
      vecteur1[j]=runif(0,1);
    }
    for(iteration=0;iteration<nbIterations;iteration++) {
      multiplicationMatriceVecteur(nbIndividus,A,vecteur1,vecteur2);
      if(iteration!=(nbIterations-1)) {
	for(i=0;i<nbDimensions;i++) {
	  vecteur1[i]=vecteur2[i];
	}
      }
    }
    valeursPropres[dimension]=0;
    denominateur=0;
    for(i=0;i<nbDimensions;i++) {
      valeursPropres[dimension]+=vecteur1[i]*vecteur2[i];
      denominateur+=vecteur1[i]*vecteur1[i];
    }
    valeursPropres[dimension]/=denominateur;
    
    /*on cherche le vecteur propre associe*/
    calculVecteurPropre(nbIndividus,A,valeursPropres[dimension],vecteursPropres[dimension]);
  
    /*for(i=0;i<nbIndividus;i++) { 
      for(j=0;j<nbIndividus;j++) {
	printf("%f ",A[i][j]);
      }
      printf("\n");
    }
    */    

    /*for(i=0;i<nbIndividus;i++) {
      printf("vp : %f\n",vecteursPropres[dimension][i]);
      }*/

    /*deflation de Wielandt*/
    /*on cherche premiere coordonnee non nulle du vecteur propre*/
    int indice;
    double norme;
    norme=0;
    for(i=0;i<nbIndividus;i++) {
      norme+=vecteursPropres[dimension][i]*vecteursPropres[dimension][i];
      if(fabs(vecteursPropres[dimension][i])>ZERO_LIMIT) {
	indice=i;
      }
    }
    norme=sqrt(norme);

    for(i=0;i<nbIndividus;i++) {
      vecteur1[i]=0;
      vecteursPropres[dimension][i]/=norme;
    }
    vecteur1[indice]=valeursPropres[dimension]/vecteursPropres[dimension][indice];

    for(i=0;i<nbIndividus;i++) {
      for(j=0;j<nbIndividus;j++) {
	A[i][j]-=vecteursPropres[dimension][i]*vecteur1[j];
      }
    }
  }

  for(dimension=1;dimension<nbDimensionsSpectrales+1;dimension++) {
        valeursPropres[dimension-1]=1.0/valeursPropres[dimension]-1.0;
    printf("val prop : %f\n",valeursPropres[dimension-1]);
    calculVecteurPropre(nbIndividus,laplacien,valeursPropres[dimension-1],vecteursPropres[dimension-1]);
  }
   
  /*projection sur nbDimensionsSpectral premiers vecteurs propres sans compter le premier*/
  for(i=0;i<nbIndividus;i++) {
    for(j=0;j<nbDimensionsSpectrales;j++) {
      individus[i]->valeursTraitees[j]=vecteursPropres[j][i];
    }
  }
  
  int nbIterationsMax,nbClustersMax;
  nbClustersMax=-1;
  nbIterationsMax=60;
  
  /*clustering par mixture models*/
  mixtureModel_t mixtureModel;
  mixtureModel.nbClusters=*nbClusters;
  mixtureModel.destroy=0;
  clusteringMixtureModels(nbIndividus,nbDimensionsSpectrales,individus,
			  &mixtureModel,nbIterationsMax,nbClustersSelection,
			  COORDONNEES,diag,nbClustersMax,NON);

  *nbClusters=mixtureModel.nbClusters;

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++) {
    ckfree(laplacien[i]);
    ckfree(A[i]);
    ckfree(B[i]);
  }
  for(i=0;i<nbDimensionsSpectrales+1;i++) {
    ckfree(vecteursPropres[i]);
  }
  ckfree(laplacien);
  ckfree(vecteursPropres);
  ckfree(valeursPropres);
  ckfree(A);
  ckfree(B);
  ckfree(vecteur1);
  ckfree(vecteur2);
  /*fin desallocation memoire*/
}
