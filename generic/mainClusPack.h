#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <unistd.h>
#include <sys/types.h>

#define MATHLIB_STANDALONE 1

#include <time.h>
#include "Rmath.h"
#include "tcl.h"
#include <limits.h>

/* Compatibilite Tcl 8.6 / Tcl 9 : les longueurs (listes, chaines) sont des
   Tcl_Size en Tcl 9 (ptrdiff_t) et des int en Tcl 8. TEA >= 3.13 passe
   deja -DTcl_Size=int quand on compile avec Tcl 8. */
#ifndef TCL_SIZE_MAX
#define Tcl_GetSizeIntFromObj Tcl_GetIntFromObj
#define TCL_SIZE_MAX INT_MAX
#ifndef Tcl_Size
typedef int Tcl_Size;
#endif
#define TCL_SIZE_MODIFIER ""
#endif

/* La librairie est chargee dans une application Tcl : les printf de
   debug du code d'origine ne doivent pas polluer stdout.
   Compiler avec -DCLUSPACK_VERBOSE pour les retrouver. */
#ifndef CLUSPACK_VERBOSE
#define printf(...) ((void)0)
#endif

/* Pour le debug memoire (AddressSanitizer, valgrind) : allocateur systeme
   au lieu de l'allocateur Tcl, qui regroupe les blocs et masque les
   debordements. Compiler avec -DCLUSPACK_SYSTEM_MALLOC. */
#ifdef CLUSPACK_SYSTEM_MALLOC
#undef ckalloc
#undef ckfree
#undef ckrealloc
#define ckalloc(n) ((char *)malloc(n))
#define ckfree(p) free((char *)(p))
#define ckrealloc(p,n) ((char *)realloc((char *)(p),(n)))
#endif

/*constantes booleennes*/
#define OUI 1
#define NON 0
#define NOT_DEFINED -1
#define OK 1
#define PAS_OK 0

/*reussite d'une fonction*/
#define SUCCES 1
#define ECHEC 0
 
/*constantes du type de cluster*/
#define NORMAL 0
#define BRUIT 1
#define VIDE 2

/*constantes des types de donnees*/
#define COORDONNEES 0
#define GAUSSIAN 1
#define SPHERICAL_GAUSSIAN 2
#define SIMPLEX 3
#define ALIGNEMENT 4
#define DISTANCES 5
#define SIMILARITES 6
#define ARBRE 7

/*constantes des programmes de clustering*/
#define MIXTURE_MODELS 0
#define BIONJ 1
#define WARD 2
#define KMEANS 3
#define KERNEL_KMEANS 4
#define KERNEL_MIXTURE_MODELS 5
#define SPECTRAL 6
#define EIGEN_MIXTURE_MODELS 7

/*constantes des criteres de selection du nombre de groupe*/
#define DPC 0
#define SECATOR 1
#define GRAPHPC 2
#define AIC 3
#define BIC 4
#define GAP 5
#define FIXED 6

/*constantes de DPC*/
#define DENSITE1 0
#define DENSITE2 1

/*constantes pour les methodes d'estimation de parametres*/
#define MOMENTS 0
#define MAXIMUM_LIKELIHOOD 1

/*constantes de types de fichiers*/
#define FICHIER_TFA 0
#define FICHIER_MSF 1
#define FICHIER_INCONNU 2

/*constantes liees a la gestion de la memoire*/
#define TAILLE_NOM 200
#define TAILLE_CHEMIN 4096
#define TAILLE_POUBELLE 40000
#define TAILLE_MAX_LIGNE 40000

/*constantes de calcul*/
#define ZERO_LIMIT 0.00000001
#define VALEUR_ENORME DBL_MAX

/*constantes diverses*/
#define AUCUN_GROUPE -1

/* Different kinds of "N(0,1)" generators :*/
typedef enum {
    BUGGY_KINDERMAN_RAMAGE,
    AHRENS_DIETER,
    BOX_MULLER,
    USER_NORM,
    INVERSION,
    KINDERMAN_RAMAGE
} N01type;


#define ML_POSINF	DBL_MAX
#define ML_NEGINF	(-DBL_MAX)
#define ML_NAN		(-DBL_MAX*(1-1e-15))
#define NA_REAL		(-DBL_MAX)

#define ML_UNDERFLOW	0

#define TRUE 1
#define FALSE 0

#define ME_NONE		0	/*	no error */
#define ME_DOMAIN	1	/*	argument out of domain */
#define ME_RANGE	2	/*	value out of range */
#define ME_NOCONV	4	/*	process did not converge */
#define ME_PRECISION	8	/*	does not have "full" precision */
#define ME_UNDERFLOW	16	/*	and underflow occured (important for IEEE)*/

#define ML_ERROR(x)	ml_error(x)
#define ML_VALID(x)	(errno == 0)
#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

#define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x)		printf(fmt,x)

#define R_FINITE(x)    R_finite(x)


#define D__0(log_p)	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define D__1(log_p)	(log_p ? 0. : 1.)			/* 1 */
#define DT_0(lower_tail,log_p)	(lower_tail ? D__0(log_p) : D__1(log_p))	/* 0 */
#define DT_1(lower_tail,log_p)	(lower_tail ? D__1(log_p) : D__0(log_p))	/* 1 */


typedef struct individu_t {
  int id;
  char *nom;
  double *valeursBrutes;
  double *valeursTraitees;
  double *inverseBandwidth;
  int cluster;
  double *vrais;
  double vraisemblance;
  char *description;
  int marque; /*0 si pas marque et 1 si marque*/
  int nbIndividusSimilaires;
  struct individu_t **individusSimilaires;
  double valeurTemp;
  double determinant;
  struct individu_t **kNearestNeighbours;
  double *distancesPivots;
  struct noeud_t *noeud;}individu_t;

/*type d'un noeud de l'arbre phylogenetique*/
typedef struct noeud_t{
  int numero;
  int groupe;
  int qualite;
  int marquage;
  char etiquette[TAILLE_NOM];
  double poids;
  double dissimilarite;
  double distance;
  double distances[3]; /*en particulier, distances[0] est la distance au noeud superieur*/
  /*                     quand l'arbre est roote*/
  struct noeud_t *copain1;
  struct noeud_t *copain2;
  struct noeud_t *copains[3]; /*en particulier, copains[0] est le pere du noeud*/
  /*                            quand l'arbre est roote*/
  int secable;
  individu_t *individu;
}noeud_t;

typedef struct mixtureModel_t{
  int nbClusters;
  int destroy;
  struct cluster_t *clusters;
}mixtureModel_t;

typedef struct cluster_t {
  int nbIndividus;
  individu_t **individus;
  int insecable;
  struct cluster_t *precedent;
  struct cluster_t *suivant;
  struct cluster_t *fils1;
  struct cluster_t *fils2;
  struct cluster_t *pere;
  struct mixtureModel_t *mixtureModel;
  double *mu;
  double **sigma;
  double *sigmaDiag;
  double logDeterminant;
  double *parametres1D;
  double *alphas;  /*pour les distributions de Dirichlet et kernel mixture models*/
  double tau; 
  double logDensiteLoiUniforme;
  double constanteNormalisation;
  int qualite;
  double sommeKxixj;
  double *sommesDiagKxixj;
  double **sommesKxixj;
  int kNearestNeighbours; /*nombre de plus proches voisin pour le cluster*/
  individu_t **pivots;
}cluster_t;

typedef struct couple {
  int i;
  int j;}couple_t;


/*****************************/
/*                           */
/*Procedure de message d'aide*/
/*                           */
/*****************************/

void helpMessageClusPack();
void DesallocArbre(struct noeud_t *arbre);
void DesallocFeuille(struct noeud_t *arbre);










