#include "mainClusPack.h"
//#include "lecfi.h"
#include "filtering.h"
#include "kmeans.h"
#include "ecrfi.h"
#include "divide.h"
#include "spectralClustering.h"
#include "kernelMixtureModels.h"
#include "lecfi.h"
#include "kernelKmeans.h"
#include "treebionj.h"
#include "alignment.h"
#include "secatorbionj.h"
#include "treeward.h"
#include "treecut.h"
#include "treedpc.h"
#include "pretraitements.h"
#include "mds.h"
#include "mixtureModels.h"
#include "GAPkmeans.h"
#include "kernels.h"
#include "densityEstimation.h"
#include "eigenMixtureModels.h"
#include "tcl.h"
#include "TclDataIn.h"
#include "TclDataOut.h"

/* try to make a thread-safe extension */
//TCL_DECLARE_MUTEX(myMutex)


/**********************************************/
/*                                            */
/*Programme principal du package de clustering*/
/*                                            */
/**********************************************/
void LibereIndividus(int nbIndividus, individu_t **individus);

/* message d'erreur et retour TCL_ERROR */
static int ErreurCluspack(Tcl_Interp *interp, const char *format, const char *arg) {
  Tcl_SetObjResult(interp, Tcl_ObjPrintf(format, arg));
  return TCL_ERROR;
}

/* Tcl_ListObjGetElements avec une longueur int : en Tcl 9 les longueurs
   sont des Tcl_Size (64 bits), le code de ClusPack travaille en int. */
static int ElementsListe(Tcl_Interp *interp, Tcl_Obj *listPtr, int *nbPtr, Tcl_Obj ***elemsPtr) {
  Tcl_Size nb;
  if (Tcl_ListObjGetElements(interp,listPtr,&nb,elemsPtr) != TCL_OK) {
    return TCL_ERROR;
  }
  if (nb > INT_MAX) {
    return ErreurCluspack(interp, "%s", "data: list too long");
  }
  *nbPtr = (int)nb;
  return TCL_OK;
}

int TclcluspackObjCmd (ClientData cdata, 
		 Tcl_Interp *interp, 
		 int objc, 
		 Tcl_Obj *const objv[]) {

  /*declaration de variables*/
  int i,j,k,typeDonnees,typeDonneesInitial,clusteringMethod,nbClustersSelection,nbClustersSelected=0;
  int typeDensite,weighting,donneesStandardisees,standardisation;
  int writeCoordinates,filtrageIndividus,longueurAlignement,nbIndividus;
  int nbDimensionsOrigine,nbDimensionsDemandees,nbDimensionsUtilisees;
  int nbClusters,nbIterationsMax,nbVoisins,individusTousIdentiques,indiceIndividu;
  int bandwidthType;
  int nbDimensionsMax=20,diag,choixDistribution,methodeEstimationParametres;
  char *buf, **sequences;
  char *arbreTexte;
  double distanceMin=0.0,seuilDissimilarite=0.0,**distances=NULL,bandwidth=1.0;
  individu_t **individus=NULL;
  noeud_t *arbre=NULL;
  /* Pour Tcl */
  int index;

  static const char *const options[] = {
    "-dt", "-cm", "-nbc", "-dt1", "-dt2", "-weighting", "-standardized_data", "-standardization", "-wc", "-nbsim", "-neighb", "-dim", "-diag", "-fd", "-distri", "-estim", "-band", NULL};

  enum options {
    CPK_DATA_TYPE,CPK_CMETHOD,CPK_NCLUSTERS,CPK_DIST1,CPK_DIST2,CPK_WEIGHTING,CPK_STANDARDIZED_DATA,CPK_STANDARDIZATION,CPK_WRITE_COORDS,CPK_NBSIM,CPK_NEIGHB,CPK_DIM,CPK_DIAG,CPK_FILTER_DIST,CPK_DISTRI,CPK_ESTIM,CPK_BAND
  };
  /* options suivies d'une valeur */
  static const int avecValeur[] = {1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1};
  /*fin declaration de variables*/

  (void)bandwidthType;

  /************************************/
  /*Traitement des parametres d'entree*/ 
  /************************************/

  if (objc < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, "data ?-dt coordinates|alignment|distances|newick|gaussian|sgaussian|simplex? -cm kmeans|hierar|bionj|mm|kkmeans|kmm|emm ?-nbc dpc|secator|aic|bic|gap|number? ?option ...?");
    return TCL_ERROR;
  }

  /*initialisation de la generation aleatoire de nombres*/
  set_seed(getpid(),(int)floor(10000.0*(double)rand()/(double)RAND_MAX));

  /*parametres par defaut*/
  typeDonnees=COORDONNEES;
  weighting=NON;
  donneesStandardisees=NON;
  standardisation=NON;
  writeCoordinates=NON;
  filtrageIndividus=NON;
  typeDensite=DENSITE1;
  nbIterationsMax=-1;
  clusteringMethod=-1;
  nbVoisins=-1;
  nbDimensionsDemandees=-1;
  nbDimensionsOrigine=0;
  nbDimensionsUtilisees=0;
  bandwidthType=ABMM;
  choixDistribution=EXPONENTIAL;
  methodeEstimationParametres=MOMENTS;
  diag=NOT_DEFINED;
  nbClusters=0;
  nbClustersSelection=DPC;
  nbIndividus=0;
  /*fin parametres par defaut*/

  i=2;
  while (i<objc) {      
    if (Tcl_GetIndexFromObj(interp, objv[i], options, "option", 0, &index) != TCL_OK) {
      return TCL_ERROR;
    }
    if (avecValeur[index]) {
      if (i+1 >= objc) {
	return ErreurCluspack(interp, "missing value for option \"%s\"", options[index]);
      }
      i++;
      buf=Tcl_GetString(objv[i]);
    } else {
      buf=NULL;
    }
    switch ((enum options) index) {
    case CPK_DATA_TYPE:
      if(strcmp(buf,"coordinates")==0) {
	typeDonnees=COORDONNEES;
      } else if(strcmp(buf,"alignment")==0) {
	typeDonnees=ALIGNEMENT; 
      } else if(strcmp(buf,"distances")==0) {
	typeDonnees=DISTANCES;
      } else if(strcmp(buf,"newick")==0) {
	typeDonnees=ARBRE;
      } else if(strcmp(buf,"gaussian")==0) {
	typeDonnees=GAUSSIAN;
      } else if(strcmp(buf,"sgaussian")==0) {
	typeDonnees=SPHERICAL_GAUSSIAN;
      } else if(strcmp(buf,"simplex")==0) {
	typeDonnees=SIMPLEX;
      } else {
	return ErreurCluspack(interp, "bad data type \"%s\": must be coordinates, alignment, distances, newick, gaussian, sgaussian or simplex", buf);
      }
      break;
    case CPK_CMETHOD:
      if(strcmp(buf,"mm")==0) {
	clusteringMethod=MIXTURE_MODELS;
      } else if(strcmp(buf,"kmeans")==0) {
	clusteringMethod=KMEANS;
      } else if(strcmp(buf,"hierar")==0) {
	clusteringMethod=WARD;
      } else if(strcmp(buf,"bionj")==0) {
	clusteringMethod=BIONJ;
      } else if(strcmp(buf,"kkmeans")==0) {
	clusteringMethod=KERNEL_KMEANS;
      } else if(strcmp(buf,"kmm")==0) {
	clusteringMethod=KERNEL_MIXTURE_MODELS;
      } else if(strcmp(buf,"emm")==0) {
	clusteringMethod=EIGEN_MIXTURE_MODELS;
      } else {
	return ErreurCluspack(interp, "bad clustering method \"%s\": must be kmeans, hierar, bionj, mm, kkmeans, kmm or emm", buf);
      }
      break;
    case CPK_NCLUSTERS:
      if(strcmp(buf,"dpc")==0) {
	nbClustersSelection=DPC;
      } else if(strcmp(buf,"secator")==0) {
	nbClustersSelection=SECATOR;
      } else if(strcmp(buf,"aic")==0) {
	nbClustersSelection=AIC;
      } else if(strcmp(buf,"bic")==0) {
	nbClustersSelection=BIC;
      } else if(strcmp(buf,"gap")==0) {
	nbClustersSelection=GAP;
      } else {
	if (Tcl_GetIntFromObj(NULL,objv[i],&nbClustersSelected) != TCL_OK || nbClustersSelected < 1) {
	  return ErreurCluspack(interp, "bad number of clusters \"%s\": must be dpc, secator, aic, bic, gap or a positive integer", buf);
	}
	nbClusters=nbClustersSelected;
	nbClustersSelection=FIXED;
      }
      break;
    case CPK_BAND:
      /* option acceptee mais ignoree (bandwidth adaptative ABMM) */
      break;
    case CPK_DIST1:
      typeDensite=DENSITE1;
      break;
    case CPK_DIST2:
      typeDensite=DENSITE2;
      break;
    case CPK_WEIGHTING:
      weighting=OUI;
      break;
    case CPK_STANDARDIZED_DATA:
      donneesStandardisees=OUI;
      break;
    case CPK_STANDARDIZATION:
      standardisation=OUI;
      donneesStandardisees=OUI;
      break;
    case CPK_WRITE_COORDS:
      writeCoordinates=OUI;
      break;
    case CPK_NBSIM:
      if (Tcl_GetIntFromObj(interp,objv[i],&nbIterationsMax) != TCL_OK) {
	return TCL_ERROR;
      }
      break;
    case CPK_NEIGHB:
      if (Tcl_GetIntFromObj(interp,objv[i],&nbVoisins) != TCL_OK) {
	return TCL_ERROR;
      }
      break;
    case CPK_DIM:
      if (Tcl_GetIntFromObj(interp,objv[i],&nbDimensionsDemandees) != TCL_OK) {
	return TCL_ERROR;
      }
      if (nbDimensionsDemandees < 1) {
	return ErreurCluspack(interp, "bad dimension \"%s\": must be a positive integer", buf);
      }
      break;
    case CPK_DIAG:
      if (Tcl_GetBooleanFromObj(interp,objv[i],&diag) != TCL_OK) {
	return TCL_ERROR;
      }
      diag = diag ? OUI : NON;
      break;
    case CPK_FILTER_DIST:
      filtrageIndividus=OUI;
      if (Tcl_GetDoubleFromObj(interp,objv[i],&distanceMin) != TCL_OK) {
	return TCL_ERROR;
      }
      break;
    case CPK_DISTRI:
      if (strcmp(buf,"e")==0) {
	choixDistribution=EXPONENTIAL;
      } else if (strcmp(buf,"g")==0){
	choixDistribution=GAMMA;
      } else if (strcmp(buf,"b")==0){
	choixDistribution=BETA;
      } else {
	return ErreurCluspack(interp, "bad distribution \"%s\": must be e, g or b", buf);
      }
      break;
    case CPK_ESTIM:
      if(strcmp(buf,"moments")==0) {
	methodeEstimationParametres=MOMENTS;
      } else if (strcmp(buf,"maxlike")==0){
	methodeEstimationParametres=MAXIMUM_LIKELIHOOD;
      } else {
	return ErreurCluspack(interp, "bad estimation method \"%s\": must be moments or maxlike", buf);
      }
      break;
    }
    i++;
  }

  /* coherence des options */
  if (clusteringMethod==-1 && typeDonnees!=ARBRE) {
    return ErreurCluspack(interp, "%s", "option -cm is required");
  }
  if (clusteringMethod==KMEANS && nbClustersSelection!=FIXED
      && nbClustersSelection!=DPC && nbClustersSelection!=GAP) {
    return ErreurCluspack(interp, "%s", "-cm kmeans requires -nbc dpc, gap or a number");
  }
  if (clusteringMethod==WARD && nbClustersSelection!=FIXED
      && nbClustersSelection!=DPC && nbClustersSelection!=SECATOR) {
    return ErreurCluspack(interp, "%s", "-cm hierar requires -nbc dpc, secator or a number");
  }
  if (clusteringMethod==MIXTURE_MODELS && nbClustersSelection!=FIXED
      && nbClustersSelection!=AIC && nbClustersSelection!=BIC) {
    return ErreurCluspack(interp, "%s", "-cm mm requires -nbc aic, bic or a number");
  }
  if (clusteringMethod==BIONJ && nbClustersSelection!=FIXED
      && nbClustersSelection!=SECATOR) {
    return ErreurCluspack(interp, "%s", "-cm bionj requires -nbc secator or a number");
  }
  if ((clusteringMethod==KERNEL_KMEANS || clusteringMethod==EIGEN_MIXTURE_MODELS)
      && nbClustersSelection!=FIXED) {
    return ErreurCluspack(interp, "%s", "-cm kkmeans and -cm emm require -nbc <number>");
  }
  if (clusteringMethod==BIONJ && typeDonnees!=ALIGNEMENT && typeDonnees!=DISTANCES) {
    return ErreurCluspack(interp, "%s", "-cm bionj requires -dt alignment or -dt distances");
  }
  typeDonneesInitial=typeDonnees;
  
  /****************************************/
  /* Fin traitement des parametres d'entree*/
  /****************************************/
  
  /*********************/
  /*                   */
  /* LECTURE DES DONNEES*/
  /*                   */
  /*********************/
  
  /*lecture des donnees*/
  if(typeDonnees==ALIGNEMENT) {
    Tcl_Obj **OInd;
    
    if (ElementsListe(interp,objv[1],&nbIndividus,&OInd) != TCL_OK) {
      return TCL_ERROR;
    }
    if (nbIndividus%2 != 0 || nbIndividus < 4) {
      return ErreurCluspack(interp, "%s", "alignment: data must be a list {name1 seq1 name2 seq2 ...} with at least 2 sequences");
    }
    nbIndividus /= 2;
    
    if (TraitementAlignement(interp,
			     OInd,
			     nbIndividus,
			     &longueurAlignement,
			     &individus,
			     &sequences) != TCL_OK) {
      return TCL_ERROR;
    }
    
    if (clusteringMethod==BIONJ) {
      /*bionj travaille directement sur les distances entre les sequences*/
      nbDimensionsUtilisees = nbIndividus;
      calculDistancesSequences(nbIndividus,
			       individus,
			       longueurAlignement,
			       sequences);
    } else {
      if (nbDimensionsDemandees==-1) {
	nbDimensionsUtilisees=2;
      } else { 
	nbDimensionsUtilisees=nbDimensionsDemandees;
      }
      if (nbVoisins==-1) {
	nbVoisins=nbDimensionsUtilisees+1;
      }
      
      /*calcul des pourcentages d'identite entre les sequences*/
      
      /*allocation memoire*/
      distances=(double **)ckalloc(sizeof(double *)*nbIndividus);
      for (i=0;i<nbIndividus;i++) {
	distances[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
      }
      /*fin allocation memoire*/
      calculDistancesApartSequences(nbIndividus,
				    longueurAlignement,
				    sequences,
				    distances);
      
      distancesToCoordonnees(nbIndividus,
			     individus,
			     distances,
			     nbDimensionsUtilisees);
      
      /*  for(i=0;i<nbIndividus;i++) {
	  for(j=0;j<nbDimensionsUtilisees;j++) {
	  printf("%.3f ",individus[i]->valeursTraitees[j]);
	  }
	  printf("\n");
	  }
	  getchar();
      */
      /*desallocation memoire*/
      for(i=0;i<nbIndividus;i++) {
	ckfree((char *)distances[i]);
      }
      ckfree((char *)distances);
      distances=NULL;
      /*fin desallocation memoire*/
    }
    for(i=0;i<nbIndividus+1;i++) {
      ckfree((char *)sequences[i]);
    }
    ckfree((char *)sequences);
  } else if (typeDonnees==ARBRE) {
    char *textBrut;
    textBrut = Tcl_GetString(objv[1]);
    /* le parseur suppose un arbre bionj : longueurs de branches et sequence fictive */
    if (strchr(textBrut,'(')==NULL || strchr(textBrut,':')==NULL) {
      return ErreurCluspack(interp, "%s", "newick: data must be a tree with branch lengths, e.g. ((a:0.1,b:0.2):0.3,c:0.1,fictive_sequence:1.0);");
    }
    if (strstr(textBrut,"fictive_sequence")==NULL) {
      return ErreurCluspack(interp, "%s", "newick: the tree must contain a leaf named \"fictive_sequence\" (outgroup used for rooting)");
    }
    arbreTexte = (char *)ckalloc(strlen(textBrut)*sizeof(char)+1);
    strcpy(arbreTexte,textBrut);
    
    /*lecture de l'arbre phylogenetique*/
    //lulu
    //lectureArbrePhylogenetique(fichierEntree,&arbreTexte);
    
    /*(re)construction de l'arbre phylogenetique*/
    
    constructionArbre(&nbIndividus,
		      arbreTexte,
		      &arbre);
    /* arbreTexte est libere par constructionArbre */
    
    nbIndividus--;
    nbDimensionsUtilisees=nbIndividus;
    individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
    for (i=0;i<nbIndividus;i++) {
      individus[i]=(individu_t *)ckalloc(sizeof(individu_t));
      individus[i]->valeursTraitees=(double *)ckalloc(sizeof(double)*nbIndividus);
      individus[i]->valeursBrutes = (double *)ckalloc(sizeof(double)*nbIndividus);
      individus[i]->cluster=-1;
      individus[i]->nom = NULL; /* rempli par associationIndividusNoeuds2 */
      individus[i]->description=NULL;
      individus[i]->id=i;
      individus[i]->cluster=AUCUN_GROUPE;
      individus[i]->nbIndividusSimilaires=0;
    }
    
    /*enracinement de l'arbre*/ 
    if (enracinementArbre(nbIndividus,&arbre) != SUCCES) {
      LibereIndividus(nbIndividus,individus);
      return ErreurCluspack(interp, "%s", "newick: the tree must contain a leaf named \"fictive_sequence\" (outgroup used for rooting)");
    }
    
    /*association des individus aux noeuds*/
    indiceIndividu=0;
    associationIndividusNoeuds2(&indiceIndividu,
				individus,
				arbre);
    
    /*calcul des distances*/
    calculDistances(nbIndividus,individus,arbre);
    
    /*the number of clusters is found by Secator*/
    secatorBionj(nbIndividus,
		 arbre,
		 &nbClusters,
		 weighting,
		 nbClustersSelection,
		 nbClustersSelected);
  } else if (typeDonnees==DISTANCES) {
    if (nbDimensionsDemandees==-1) {
      nbDimensionsUtilisees=2;
    } else {
      nbDimensionsUtilisees=nbDimensionsDemandees;
    }
    Tcl_Obj **OInd;
    if (ElementsListe(interp,objv[1],&nbIndividus,&OInd) != TCL_OK) {
      return TCL_ERROR;
    }
    if (nbIndividus < 2) {
      return ErreurCluspack(interp, "%s", "distances: data must be a list of at least 2 rows {name d1 ... dn}");
    }
    if (TraitementDistances(interp,OInd,nbIndividus,nbDimensionsUtilisees,&individus,&distances) != TCL_OK) {
      return TCL_ERROR;
    }
    /*les distances sont recopiees dans les individus (filtrage, bionj) ;
      la matrice est liberee, et reconstruite apres filtrage si besoin*/
    for(i=0;i<nbIndividus;i++) {
      for(j=0;j<nbIndividus;j++) {
	individus[i]->valeursTraitees[j]=distances[i][j];
      }
      ckfree((char *)distances[i]);
    }
    ckfree((char *)distances);
    distances=NULL;
    if (clusteringMethod==BIONJ) {
      nbDimensionsUtilisees=nbIndividus;
    }
  } else {
    /* COORDONNEES - simplex */

    Tcl_Size nbCol;
    Tcl_Obj **OInd;

    if (ElementsListe(interp,objv[1],&nbIndividus,&OInd) != TCL_OK) {
      return TCL_ERROR;
    }
    if (nbIndividus < 2) {
      return ErreurCluspack(interp, "%s", "coordinates: data must be a list of at least 2 rows {name x1 x2 ...}");
    }
    if (Tcl_ListObjLength(interp,OInd[0],&nbCol) != TCL_OK) {
      return TCL_ERROR;
    }
    if (nbCol < 2) {
      return ErreurCluspack(interp, "%s", "coordinates: each row must be {name x1 x2 ...}");
    }
    individus = (individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
    nbDimensionsOrigine = (int)(nbCol-1);    
    for(i=0;i<nbIndividus;i++) {
      individus[i]=(individu_t *)ckalloc(sizeof(individu_t));
      individus[i]->id=i;
      individus[i]->cluster=AUCUN_GROUPE;
      individus[i]->nbIndividusSimilaires=0;
      individus[i]->cluster=-1;
      individus[i]->valeursBrutes = (double *)ckalloc(sizeof(double)*nbDimensionsOrigine);
      individus[i]->valeursTraitees = (double *)ckalloc(sizeof(double)*nbDimensionsOrigine);
      individus[i]->nom=NULL;
      individus[i]->description=NULL;
    }

    if (TraitementCoordonnees(interp,
			      OInd,
			      nbIndividus,
			      nbDimensionsOrigine,
			      &individus,
			      typeDonnees) != TCL_OK) {
      LibereIndividus(nbIndividus,individus);
      return TCL_ERROR;
    }
    
    filtreDimensions(nbIndividus,
		     nbDimensionsOrigine,
		     &nbDimensionsUtilisees,
		     individus);
    //printf("nbDimensionsUtilisees : %d\n",nbDimensionsUtilisees);
  }
  
  /*****************************/
  /*                           */
  /* PRE-TRAITEMENT DES DONNEES*/
  /*                           */
  /*****************************/
  
  if (filtrageIndividus==OUI) {
    //fprintf(stderr,"filtrage\n");
    //fflush(stderr);
    /*filtrage des individus suivant la valeur seuil donnee et le type des donnees*/
    filtreIndividus(&nbIndividus,nbDimensionsUtilisees,individus,typeDonnees,distanceMin,
		    (typeDonnees==DISTANCES || clusteringMethod==BIONJ) ? OUI : NON);
  }
  if((typeDonnees==COORDONNEES)&&(nbDimensionsDemandees>0)) {
    /*reduction de dimension*/
    reductionDimension(nbIndividus,individus,nbDimensionsOrigine,nbDimensionsDemandees);
    nbDimensionsUtilisees=nbDimensionsDemandees;
  } else if(typeDonnees==DISTANCES && clusteringMethod!=BIONJ) {
    /*conversion des distances en coordonnees*/
    if(nbDimensionsDemandees==-1) {
      nbDimensionsUtilisees=2;
    } else {
      nbDimensionsUtilisees=nbDimensionsDemandees;
    }
    /*matrice des distances entre les individus restant apres filtrage*/
    distances=(double **)ckalloc(sizeof(double *)*nbIndividus);
    for(i=0;i<nbIndividus;i++) {
      distances[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
      for(j=0;j<nbIndividus;j++) {
	distances[i][j]=individus[i]->valeursTraitees[individus[j]->id];
      }
    }
    distancesToCoordonnees(nbIndividus,individus,distances,nbDimensionsUtilisees);
    for(i=0;i<nbIndividus;i++) {
      ckfree((char *)distances[i]);
    }
    ckfree((char *)distances);
    distances=NULL;
  }
  
  if (standardisation==OUI) {
    //fprintf(stderr,"standardisation\n");
    //fflush(stderr);
    /*on normalise les individus*/   
    i=standardisationIndividus(nbIndividus,nbDimensionsUtilisees,individus);
    if (i>=0) {
      Tcl_SetObjResult(interp,Tcl_ObjPrintf("-standardization: standard deviation of \"%s\" is zero",individus[i]->nom));
      LibereIndividus(nbIndividus,individus);
      return TCL_ERROR;
    }
    if(clusteringMethod==MIXTURE_MODELS) {
      //printf("donnees spheriques\n");
      typeDonnees=SPHERICAL_GAUSSIAN;
      projectionSurSphere(nbIndividus,nbDimensionsUtilisees,individus);
      nbDimensionsUtilisees--;	
    }
  }
  
  /***************/
  /*             */
  /* CLUSTERING  */
  /*             */
  /***************/
  if(diag==NOT_DEFINED) {
    /*covariance pleine : chaque cluster doit contenir au moins d+1 individus
      (kmixtureModels), impossible si nbIndividus <= d, meme pour k=1
      => covariance diagonale*/
    if((nbDimensionsUtilisees>nbDimensionsMax)||(nbIndividus<nbDimensionsUtilisees+1)) {
      diag=OUI;
    } else {
      diag=NON;
    }
  }
  
  /*clustering method : KMEANS */
  if (clusteringMethod==KMEANS) {
    /*the number of clusters is fixed*/
    if (nbClustersSelection==FIXED) {
      kmeans(nbIndividus,nbDimensionsUtilisees,individus,nbClustersSelected);
      nbClusters=nbClustersSelected;
    } else if ((nbClustersSelection==DPC)||(nbClustersSelection==GRAPHPC)) {
      /*the number of clusters is automatically found by DPC*/
      
      divide(nbIndividus,
	     nbDimensionsUtilisees,
	     &nbClusters,
	     individus,
	     typeDonnees,
	     clusteringMethod,
	     nbClustersSelection,
	     donneesStandardisees,
	     typeDensite,
	     nbVoisins);
    } else if(nbClustersSelection==GAP) {
      GAPkmeans(nbIndividus,nbDimensionsUtilisees,individus,-1,&nbClusters);
    }
  } else if (clusteringMethod==KERNEL_KMEANS) {
    int choixNoyau=RADIAL_GAUSSIAN;
    
    /*2026-10-08 : bandwidth par calculPilotAdaptiveBandwidth (bandwidth
      global + densite pilote a noyau) au lieu de mixtureModelDensityEstimation
      + calculMMDensityBandwidth, qui donnaient des bandwidths aberrants
      (voir kernels.c) ; le modele de melange n'etait pas libere.*/
    if (nbClustersSelected>nbIndividus) {
      Tcl_SetObjResult(interp, Tcl_ObjPrintf("-cm kkmeans: %d clusters requested for %d individuals", nbClustersSelected, nbIndividus));
      LibereIndividus(nbIndividus,individus);
      return TCL_ERROR;
    }
    /*clustering method : KMEANS */
    /*allocation memoire*/
    for(i=0;i<nbIndividus;i++) {
      individus[i]->inverseBandwidth=(double *)ckalloc(sizeof(double)*nbDimensionsUtilisees);
    }
    /*fin allocation memoire*/
    
    /*adaptive bandwidth (Abramson)*/
    calculPilotAdaptiveBandwidth(nbDimensionsUtilisees,nbIndividus,individus);
    
    /*kernel-Kmeans*/
    j=kernelKmeans(nbIndividus,
		 nbDimensionsUtilisees,
		 individus,
		 nbClustersSelected,
		 choixNoyau,
		 0,
		 NULL);
    for(i=0;i<nbIndividus;i++) {
      ckfree((char *)individus[i]->inverseBandwidth);
      individus[i]->inverseBandwidth=NULL;
    }
    if (j==ECHEC) {
      Tcl_SetObjResult(interp, Tcl_NewStringObj("kernel k-means clustering failed", -1));
      LibereIndividus(nbIndividus,individus);
      return TCL_ERROR;
    }
  } else if (clusteringMethod==KERNEL_MIXTURE_MODELS) {
    int accelerationCalculs=NON;
    int choixNoyau=RADIAL_GAUSSIAN;
    int nbPivots=100;
    
    mixtureModel_t mixtureModelPourDensityEstimation;
    /*estimation de la densite par modele de melanges*/ 
    /*allocation memoire*/
    for(i=0;i<nbIndividus;i++) {
      individus[i]->inverseBandwidth=(double *)ckalloc(sizeof(double)*nbDimensionsUtilisees);
    }
    /*fin allocation memoire*/
    mixtureModelDensityEstimation(nbIndividus,nbDimensionsUtilisees,individus,0,&mixtureModelPourDensityEstimation);
    
    /*adaptive bandwidth par estimation de la densite*/
    calculMMDensityBandwidth(nbDimensionsUtilisees,nbIndividus,individus);
    
    /*the number of clusters is fixed*/
    mixtureModel_t mixtureModel;
    if(nbClustersSelection==FIXED) {
      mixtureModel.nbClusters=nbClusters;
    }
    /*nb iterations for EM*/
    if(nbIterationsMax==-1) {
      nbIterationsMax=60;
    }
    /*kernel mixture models (matrice des noyaux inutile sans acceleration)*/ 
    kernelMixtureModels(nbIndividus,
			nbDimensionsUtilisees,
			individus,
			&mixtureModel,
			nbIterationsMax,
			nbClustersSelection,
			-1,
			choixNoyau,
			choixDistribution,
			accelerationCalculs,
			nbPivots,
			methodeEstimationParametres,
			NULL);
    nbClusters=mixtureModel.nbClusters;
  } else if (clusteringMethod==EIGEN_MIXTURE_MODELS) {
    if(nbIterationsMax==-1) {
      nbIterationsMax=60;
    }
    typeDonnees=GAUSSIAN;
    eigenMixtureModels(nbIndividus,
		       nbDimensionsUtilisees,
		       individus,
		       nbIterationsMax,
		       nbClustersSelected,
		       &nbClusters,diag);
  } else if (clusteringMethod==SPECTRAL) {
    /*clustering method: SPECTRAL*/
    spectralClustering(nbIndividus,
		       nbDimensionsUtilisees,
		       individus,
		       nbClustersSelection,
		       &nbClusters,
		       bandwidthType,
		       bandwidth,
		       diag);
  } else if (clusteringMethod==WARD) {
    /*clustering method : WARD = hierarchic*/
    noeud_t **GrosArbre;
    noeud_t *TasFeuilles;
    TasFeuilles = classificationWard(nbIndividus,
		       nbDimensionsUtilisees,
		       individus,
		       typeDonnees,
		       &GrosArbre);
    arbre=GrosArbre[0];
    ckfree((char *)GrosArbre);
    /*the number of clusters is fixed*/
    if (nbClustersSelection==FIXED) {
      nbClustersToDissimilarityThreshold(nbIndividus,
					 nbClustersSelected,
					 arbre,
					 &seuilDissimilarite);
    } else if (nbClustersSelection==SECATOR) {
      /*the number of clusters is found by Secator*/
      secatorTree(nbIndividus,
		  arbre,
		  &seuilDissimilarite,
		  &nbClusters);
    }
    if ((nbClustersSelection==FIXED)||(nbClustersSelection==SECATOR)) {
      /*creation des groupes*/
      decoupageArbreSeuilDissimilarite(nbIndividus,
				       arbre,
				       seuilDissimilarite,
				       &nbClusters);
    } else {
      decoupageArbreDPC(nbIndividus,
			nbDimensionsUtilisees,
			individus,
			arbre,
			&nbClusters,
			typeDonnees,
			nbClustersSelection,
			donneesStandardisees,
			typeDensite,
			nbVoisins);
    }

    DesallocArbre(arbre->copain1);
    DesallocArbre(arbre->copain2);
    ckfree((char *)arbre);
    arbre=NULL;
    ckfree((char *)TasFeuilles);

  } else if (clusteringMethod==MIXTURE_MODELS) {
    /*clustering method : mixture models*/
    /*nb iterations for EM*/
    if (nbIterationsMax==-1) {
      nbIterationsMax=60;
    }
    
    if ((nbClustersSelection==FIXED)||(nbClustersSelection==AIC)||(nbClustersSelection==BIC)) {
      /*the number of clusters is fixed*/

      mixtureModel_t mixtureModel;
      //mixtureModel=(mixtureModel_t *)ckalloc(sizeof(mixtureModel_t));
      if (nbClustersSelection==FIXED) {
	mixtureModel.nbClusters=nbClusters;
      }
      
      /*seul le nombre de clusters est utilise : destroy=1 fait renvoyer ce nombre
	et libere les structures internes*/
      mixtureModel.destroy=1;
      if(typeDonnees==COORDONNEES||typeDonnees==ALIGNEMENT||typeDonnees==DISTANCES) {
	/*coordonnees (eventuellement issues du MDS)*/
	typeDonnees=GAUSSIAN;
      }
      nbClusters=clusteringMixtureModels(nbIndividus,
					 nbDimensionsUtilisees,
					 individus,
					 &mixtureModel,
					 nbIterationsMax,
					 nbClustersSelection,
					 typeDonnees,
					 diag,-1,NON);
      if (nbClusters<=0) {
	/*ECHEC : l'EM n'a converge pour aucun nombre de clusters*/
	LibereIndividus(nbIndividus,individus);
	return ErreurCluspack(interp, "%s", "mixture model clustering failed (EM did not converge)");
      }
    }
  } else if (clusteringMethod==BIONJ) {
    /*clustering method : bionj*/

    /*on verifie que les individus ne sont pas tous identiques*/
    individusTousIdentiques=OUI;
    for (i=0;i<nbIndividus;i++) {
      for(j=i+1;j<nbIndividus;j++) {
	if (individus[i]->valeursTraitees[j]>ZERO_LIMIT) {
	  individusTousIdentiques=NON;
	}
      }
    }
    
    if (individusTousIdentiques==OUI) {
      nbClusters=1;
      for (i=0;i<nbIndividus;i++) {
	individus[i]->cluster=0;
      }
    } else {
      if (creationArbreBionj(nbIndividus,
			     individus,
			     &arbre) != SUCCES) {
	LibereIndividus(nbIndividus,individus);
	return ErreurCluspack(interp, "%s", "bionj: unable to build the tree (temporary files)");
      }
      /*the number of clusters is found by Secator*/
      if ((nbClustersSelection==SECATOR)||(nbClustersSelection==FIXED)) {
	secatorBionj(nbIndividus,
		     arbre,
		     &nbClusters,
		     weighting,
		     nbClustersSelection,
		     nbClustersSelected);
      }
    }
  }
  
  /*************************/
  /*                       */
  /* ECRITURE DES RESULTATS*/
  /*                       */
  /*************************/
  if ((standardisation==OUI)&&(clusteringMethod==MIXTURE_MODELS)) {
    nbDimensionsUtilisees++;
  }


  if (writeCoordinates==OUI && typeDonneesInitial!=COORDONNEES
      && typeDonneesInitial!=GAUSSIAN && typeDonneesInitial!=SPHERICAL_GAUSSIAN
      && typeDonneesInitial!=SIMPLEX) {
    /*pas de coordonnees d'origine : on renvoie les valeurs traitees
      (coordonnees MDS, ou distances pour bionj/newick)*/
    for(i=0;i<nbIndividus;i++) {
      for(j=0;j<nbDimensionsUtilisees;j++) {
	individus[i]->valeursBrutes[j]=individus[i]->valeursTraitees[j];
      }
      for(k=0;k<individus[i]->nbIndividusSimilaires;k++) {
	for(j=0;j<nbDimensionsUtilisees;j++) {
	  individus[i]->individusSimilaires[k]->valeursBrutes[j]=
	    individus[i]->individusSimilaires[k]->valeursTraitees[j];
	}
      }
    }
  }

  if (nbDimensionsDemandees>0 && typeDonneesInitial==COORDONNEES) {
    RetourneClusters(interp,
		     nbIndividus,
		     nbDimensionsOrigine,
		     individus,
		     nbClusters,
		     nbClustersSelection,
		     writeCoordinates);
  } else {
    RetourneClusters(interp,
		     nbIndividus,
		     nbDimensionsUtilisees,
		     individus,
		     nbClusters,
		     nbClustersSelection,
		     writeCoordinates);
  }

  /*l'arbre bionj/newick est libere par secatorBionj*/
  LibereIndividus(nbIndividus,individus);

  return TCL_OK;
}

/* liberation des individus (y compris ceux filtres comme similaires) */
static void LibereIndividu(individu_t *ind) {
  int k;
  if (ind == NULL) return;
  for (k=0;k<ind->nbIndividusSimilaires;k++) {
    LibereIndividu(ind->individusSimilaires[k]);
  }
  if (ind->nbIndividusSimilaires>0) {
    ckfree((char *)ind->individusSimilaires);
  }
  ckfree((char *)ind->nom);
  ckfree((char *)ind->description);
  ckfree((char *)ind->valeursTraitees);
  ckfree((char *)ind->valeursBrutes);
  ckfree((char *)ind);
}

void LibereIndividus(int nbIndividus, individu_t **individus) {
  int i;
  if (individus == NULL) return;
  for(i=0;i<nbIndividus;i++) {
    LibereIndividu(individus[i]);
  }
  ckfree((char *)individus);
}

/*
void helpMessageClusPack() {
  printf("Program usage : cluspack file options\n\n");
  printf("********************OPTIONS BELOW********************\n\n");
  printf("-dt=[coordinates|gaussian|sgaussian|newick|alignment|distances] (dt stands for data_type)\n");
  printf("-cm=[kmeans|hierar|bionj|mm|spectral] (cm stands for clusteringMethod, mm for mixture models))\n");
  printf("-nbc=[secator|dpc|aic|bic|number]  (nbc stands for method for computing the number \
of cluster and number is really a number like 4 etc.)\n");
  printf("-neighb=[number of neighbours] (for estimating the density for DPC, equal to dimension+1 by default)\n");
  printf("-dim=[dimension] (for dimension reduction when dealing with alignments and using DPC equanl to 2 by default)\n");
  printf("[-dt1|-dt2][-standardization][-standardized_data][-wc] (dt1 stands for density1 and wc for writeCoordinates)\n");
  printf("[-nbsim=nbsimulations]\n");
  printf("[-otfa=outputFile for alignment]\n");
  printf("[-oclu=outputFile for clustering]\n");

  exit(0);
}
*/

int DLLEXPORT Tclcluspack_Init(Tcl_Interp *interp) {
  if (
#ifdef USE_TCL_STUBS 
      Tcl_InitStubs(interp, "8.5-10", 0)
#else
      Tcl_PkgRequire(interp, "Tcl", "8.5-10", 0)
#endif /* USE_TCL_STUBS */
      == NULL) {
    return TCL_ERROR;
  }
  
  Tcl_CreateObjCommand(interp,"Tclcluspack",TclcluspackObjCmd,NULL,NULL);
  
  if (Tcl_PkgProvide(interp,"tclcluspack","0.1") != TCL_OK) {
    return TCL_ERROR;
  }

  return TCL_OK;
}


void DesallocArbre(noeud_t *arbre) 
{
  if (arbre == NULL) return;

  if (arbre->qualite == FEUILLE) {
    //fprintf(stderr,"  ==> suis dans feuile %x\n",arbre);
    //fflush(stderr);
    //DesallocFeuille((noeud_t)arbre);
    //ckfree((char *)arbre);
    return;
  }

  //fprintf(stderr,"desalloc %x ",arbre);
  //fprintf(stderr,"cop1 %x ",arbre->copain1);
  //fprintf(stderr,"cop2 %x ",arbre->copain2);
  //fprintf(stderr,"| %d\n",arbre->qualite);
  //fflush(stderr);

  if (arbre->copain1 != NULL) {
    DesallocArbre(arbre->copain1);
    arbre->copain1=NULL;
  }
  if (arbre->copain2 != NULL) {
    DesallocArbre(arbre->copain2);
    arbre->copain2=NULL;
  }

  ckfree((char *)arbre);
  return;





}

void DesallocFeuille(noeud_t *arbre) {
  ckfree((char *)arbre);
}

int DLLEXPORT Tclcluspack_SafeInit(Tcl_Interp *interp) {
  return Tclcluspack_Init(interp);
} /* Tclcluspack_SafeInit */

    
