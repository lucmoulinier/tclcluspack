#include "mainClusPack.h"
#include "lecfi.h"
#include "ecrfi.h"
#include "apprentissageKernels.h"
#include "predictionKernels.h" 
#include "kernels.h"
#include "mainKernelBayesClassifier.h"

int main(int argc,char **argv){
  int i,tailleTrainingSet,tailleTestSet,nbDimensions,nbClusters,accelerationCalculs,diag,choixNoyau;
  int bandwidthType,choixDistribution,methodeEstimationParametres;
  double bandwidth,scale,regularisation;
  cluster_t *clusters;
  individu_t **trainingSet,**testSet;
  char *buf;
  
  set_seed(getpid(),(int)floor(10000.0*(double)rand()/(double)RAND_MAX));

  if(argc<5) {
    helpMessageKernelBayesClassifier();
  }
  
  accelerationCalculs=NON;

  choixNoyau=RADIAL_GAUSSIAN;
  //choixNoyau=EPANECHNIKOV;
  //choixNoyau=POLYNOMIAL2;

  bandwidthType=ABMM;
  choixDistribution=GAMMA;
  methodeEstimationParametres=MOMENTS;
  //methodeEstimationParametres=MAXIMUM_LIKELIHOOD;
  scale=1.0;
  int verbose=NON;
  regularisation=-1.0;
  for(i=5;i<argc;i++) {
    if((strcmp(argv[i],"-help")==0)||(strcmp(argv[i],"-h")==0)) {
      helpMessageKernelBayesClassifier();
    } else if(strstr(argv[i],"-band=")!=NULL) {
      buf=strstr(argv[i],"=")+1;
      if(strcmp(buf,"fixed")==0) {
	bandwidth=1.0;
	bandwidthType=FIXED_BANDWIDTH;
      } else if(strcmp(buf,"adaAB")==0) {
	bandwidthType=AB;
      } else if(strcmp(buf,"adaABCW")==0) {
	bandwidthType=ABCW;
      } else if(strcmp(buf,"adaABMM")==0) {
	bandwidthType=ABMM;
      } else if(strcmp(buf,"adaABNN")==0) {
	bandwidthType=ABNN;
      } else {
	bandwidth=atof(buf);
	printf("bandwidth : %f\n",bandwidth);
	bandwidthType=FIXED_BANDWIDTH;
      }
    } else if(strstr(argv[i],"-scale=")!=NULL) {
      buf=strstr(argv[i],"=")+1;
      scale=atof(buf);
    } else if(strstr(argv[i],"-kernel=")!=NULL) {
      buf=strstr(argv[i],"=")+1;
      if(strcmp(buf,"rg")==0) {
	choixNoyau=RADIAL_GAUSSIAN;
      } else if(strcmp(buf,"e")==0) {
	choixNoyau=EPANECHNIKOV;
      }else if(strcmp(buf,"poly2")==0) {
	choixNoyau=POLYNOMIAL2;
      } else if(strcmp(buf,"poly3")==0) {
	choixNoyau=POLYNOMIAL3;
      } else if(strcmp(buf,"poly4")==0) {
	choixNoyau=POLYNOMIAL4;
      } else {
	printf("kernel %s does not existing\n",buf);
	exit(0);
      }
    } else if(strcmp(argv[i],"-boost")==0) {
      accelerationCalculs=OUI;
    } else if(strstr(argv[i],"-boost=")!=NULL) {
      accelerationCalculs=OUI;
      buf=strstr(argv[i],"=")+1;
    } else if(strstr(argv[i],"-estim=")!=NULL) {
	buf=strstr(argv[i],"=")+1;
	if(strcmp(buf,"moments")==0) {
	  methodeEstimationParametres=EXPONENTIAL;
	} else if (strcmp(buf,"maxlike")==0){
	  methodeEstimationParametres=MAXIMUM_LIKELIHOOD;
	} else {
	  printf("optimization method %s does not existing\n",buf);
	  exit(0);
	}
    } else if(strstr(argv[i],"-verbose")!=NULL) {
      verbose=OUI;
    } else if(strstr(argv[i],"-reg=")!=NULL) {
	buf=strstr(argv[i],"=")+1;
	regularisation=atof(buf);
    }
  }
  
  /*lecture des coordonnees du training set*/
  lectureFichierCoordonnees(&tailleTrainingSet,&nbDimensions,&trainingSet,argv[1],COORDONNEES);
  
  /*lecture des clusters du training set*/
  lectureEtiquettes(tailleTrainingSet,trainingSet,argv[2]);

  /*lecture des coordonnees du test set*/
  lectureFichierCoordonnees(&tailleTestSet,&nbDimensions,&testSet,argv[3],COORDONNEES);

  /*apprentissage sur le training set*/
  mixtureModel_t mixtureModelPourDensite,**mixtureModels;
  apprentissageKernels(&nbDimensions,&tailleTrainingSet,trainingSet,tailleTestSet,testSet,&nbClusters,
		       &clusters,&mixtureModels,&mixtureModelPourDensite,bandwidthType,bandwidth,accelerationCalculs,
		       &diag,choixNoyau,&choixDistribution,scale,methodeEstimationParametres,verbose,regularisation);

  /*prediction sur le test set*/
  predictionKernels(nbDimensions,nbClusters,clusters,tailleTrainingSet,trainingSet,tailleTestSet,testSet,mixtureModels,
		    mixtureModelPourDensite,bandwidthType,bandwidth,accelerationCalculs,diag,choixNoyau,choixDistribution,
		    methodeEstimationParametres);

  /*affichage des clusters predits*/
  ecritureFichierEtiquettes(tailleTestSet,testSet,argv[4]);

  exit(0);
}

/*****************************/
/*                           */
/*Procedure de message d'aide*/
/*                           */
/*****************************/

void helpMessageKernelBayesClassifier(){
  printf("Usage : kernelClassification coordonneesTrainingSet etiquettesTrainingSet coordonneesTestSet \
etiquettesTestSet [-band=fixed|ada1|ada2] [-kernel=rg|poly2] [-distri=iso|diagg|cg] [-boost]\n");
  exit(0);
}
