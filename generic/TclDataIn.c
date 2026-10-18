#include "mainClusPack.h"
//#include "lecfi.h"
#include <ctype.h>
#include "tcl.h"
#include "TclDataIn.h"

/****************************************************/
/*                                                  */
/*  Remplace lecfi.c de Nicolas.                    */
/*                                                  */
/* Traite donnes en entrees:                        */
/* - arbre (au format bionj)                        */
/*                                                  */
/****************************************************/

/********************************************/
/*                                          */
/*Procedure de lecture du fichier contenant */
/*les coordonnees des individus a classifier*/
/*Les individus sont stockes dans un tableau*/
/*de pointeurs sur structure                */
/*                                          */
/********************************************/

int TraitementCoordonnees(Tcl_Interp *interp,
			  Tcl_Obj **OLinePtr,
			  int nbIndividus,
			  int nbDimensions,
			  individu_t ***individus, 
			  int typeDonnees)
{
  /*declaration des variables*/
  int i,j;
  Tcl_Size nbCol,lgt;
  char *nomIndividu;
  double somme;
  Tcl_Obj **OPtr;
  /*fin declaration des variables*/

  /*
    if (Tcl_ListObjLength(interp, OLinePtr[0],&nbCol) != TCL_OK) {
    return TCL_ERROR;
    }
    *nbDimensions = nbCol-1;
    */

  /* allocation memoire */
  //*individus = (individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);

  /*
  for(i=0;i<nbIndividus;i++) {
    (*individus)[i]=(individu_t *)ckalloc(sizeof(individu_t));
    (*individus)[i]->id=i;
    (*individus)[i]->cluster=AUCUN_GROUPE;
    (*individus)[i]->valeursBrutes = (double *)ckalloc(sizeof(double)*(*nbDimensions));
    (*individus)[i]->valeursTraitees = (double *)ckalloc(sizeof(double)*(*nbDimensions));
    (*individus)[i]->nbIndividusSimilaires=0;
    (*individus)[i]->cluster=-1;
  }
  */
  /*fin allocation memoire*/

  /* Assignation */
  somme=0.0;
  for(i=0;i<nbIndividus;i++) {
    if (Tcl_ListObjGetElements(interp, OLinePtr[i], &nbCol,&OPtr) != TCL_OK) {
      return TCL_ERROR;
    }
    if (nbCol != nbDimensions+1) {
      Tcl_SetObjResult(interp,Tcl_ObjPrintf(
        "coordinates: row %d has %" TCL_SIZE_MODIFIER "d values, expected %d",i,nbCol-1,nbDimensions));
      return TCL_ERROR;
    }

    for(j=0;j<nbCol;j++) {
      if (j==0) {
	nomIndividu = Tcl_GetStringFromObj(OPtr[j],&lgt);
	lgt+=2;
	(*individus)[i]->nom = (char *)ckalloc(sizeof(char)*lgt);
	(*individus)[i]->description = (char *)ckalloc(sizeof(char)*lgt);
	strcpy((*individus)[i]->nom,nomIndividu);
	strcpy((*individus)[i]->description,nomIndividu);

	somme=0.0;
      } else {
	if (Tcl_GetDoubleFromObj(interp, OPtr[j], &((*individus)[i]->valeursBrutes[j-1])) != TCL_OK) {
	  return TCL_ERROR;
	}
	(*individus)[i]->valeursTraitees[j-1]=(*individus)[i]->valeursBrutes[j-1];
	if (typeDonnees==SIMPLEX) {
	  if ((*individus)[i]->valeursTraitees[j-1]<ZERO_LIMIT)
	    (*individus)[i]->valeursTraitees[j-1]+=0.001;
	  somme+=(*individus)[i]->valeursTraitees[j-1];
	}
      }	  
    }
    if (typeDonnees==SIMPLEX) {
      for (j=0;j<nbDimensions;j++) {
	(*individus)[i]->valeursTraitees[j]/=somme;
      }
    }
  }

  return TCL_OK;
}


/********************************************/
/*                                          */
/* Matrice de distances : une liste Tcl de  */
/* nbIndividus lignes {nom d1 d2 ... dn}    */
/* Les valeurs des individus sont allouees  */
/* sur max(n,nbDimensions) pour pouvoir     */
/* contenir soit des coordonnees soit les   */
/* distances (bionj)                        */
/*                                          */
/********************************************/

int TraitementDistances(Tcl_Interp *interp,
			Tcl_Obj **OLinePtr,
			int nbIndividus,
			int nbDimensions,
			individu_t ***individus,
			double ***distances)
{
  int i,j,taille;
  Tcl_Size nbCol;
  char *nom;
  Tcl_Obj **OPtr;

  taille=(nbIndividus>nbDimensions)?nbIndividus:nbDimensions;

  *individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  *distances=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++) {
    (*individus)[i]=(individu_t *)ckalloc(sizeof(individu_t));
    (*individus)[i]->id=i;
    (*individus)[i]->cluster=AUCUN_GROUPE;
    (*individus)[i]->nbIndividusSimilaires=0;
    (*individus)[i]->valeursBrutes=(double *)ckalloc(sizeof(double)*taille);
    (*individus)[i]->valeursTraitees=(double *)ckalloc(sizeof(double)*taille);
    (*individus)[i]->nom=NULL;
    (*individus)[i]->description=NULL;
    (*distances)[i]=(double *)ckalloc(sizeof(double)*nbIndividus);
  }

  for(i=0;i<nbIndividus;i++) {
    if (Tcl_ListObjGetElements(interp,OLinePtr[i],&nbCol,&OPtr) != TCL_OK) {
      return TCL_ERROR;
    }
    if (nbCol != nbIndividus+1) {
      Tcl_SetObjResult(interp,Tcl_ObjPrintf(
        "distances: row %d has %" TCL_SIZE_MODIFIER "d distances, expected %d",i,nbCol-1,nbIndividus));
      return TCL_ERROR;
    }
    nom=Tcl_GetString(OPtr[0]);
    (*individus)[i]->nom=strcpy(ckalloc(strlen(nom)+1),nom);
    for(j=0;j<nbIndividus;j++) {
      if (Tcl_GetDoubleFromObj(interp,OPtr[j+1],&((*distances)[i][j])) != TCL_OK) {
	return TCL_ERROR;
      }
    }
  }
  return TCL_OK;
}


int TraitementAlignement(Tcl_Interp *interp,
			 Tcl_Obj **OLine,
			 int nbIndividus,
			 int *longueurAlignement,
			 individu_t ***individus,
			 char ***sequences)
{
  /* Declaration variables */
  int i, j, c;
  Tcl_Size lgt;
  char *seq, *nom;
  /* fin declaration des variables a*/

  /* longueur de l'alignement = longueur de la plus longue sequence */
  *longueurAlignement = 0;
  for(i=0;i<nbIndividus;i++) {
    Tcl_GetStringFromObj(OLine[2*i+1],&lgt);
    if (lgt > INT_MAX) {
      Tcl_SetObjResult(interp,Tcl_NewStringObj("alignment: sequence too long",-1));
      return TCL_ERROR;
    }
    if (lgt > *longueurAlignement) *longueurAlignement = (int)lgt;
  }

  /*allocation memoire*/
  *sequences = (char **)ckalloc(sizeof(char *)*(nbIndividus+1));
  for(i=0;i<nbIndividus+1;i++) {
    (*sequences)[i]=(char *)ckalloc(sizeof(char)*(*longueurAlignement));
  }
  *individus = (individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      (*individus)[i]=(individu_t *)ckalloc(sizeof(individu_t));
      (*individus)[i]->id = i;
      (*individus)[i]->cluster = -1;
      (*individus)[i]->valeursBrutes = (double *)ckalloc(sizeof(double)*nbIndividus);
      (*individus)[i]->valeursTraitees = (double *)ckalloc(sizeof(double)*nbIndividus);
      (*individus)[i]->description = NULL;
      (*individus)[i]->nbIndividusSimilaires = 0;
    }
  /*fin allocation memoire*/

  for(i=0;i<nbIndividus;i++) {
    nom = Tcl_GetString(OLine[2*i]);
    (*individus)[i]->nom = strcpy(ckalloc(strlen(nom)+1),nom);
    seq = Tcl_GetStringFromObj(OLine[2*i+1],&lgt);
    for(j=0;j<*longueurAlignement;j++) {
      if (j>=lgt) {
	/* sequence plus courte que la premiere : on complete */
	(*sequences)[i][j]='O';
	continue;
      }
      c = toupper((unsigned char)seq[j]);
      if ((c-'A'<0)||(c-'A'>24)) {
	(*sequences)[i][j]='O';
      } else if ((c=='B')||(c=='J')||(c=='O')||(c=='U')) {
	(*sequences)[i][j]='O';
      } else {
	(*sequences)[i][j]=c;
      }
    }
  }

  for(i=0;i<*longueurAlignement;i++)
    (*sequences)[nbIndividus][i]='P';

  return TCL_OK;
}



/*
  qu'est ce que fait dans les projets
  Ce que je fais et que les autres ne font pas.
  personnes referentes. en quoi.
  anne friedrich, post doc (on s'en fout)
  
  + structure pour l'equipe, pour l'unite, pour communaute nationale et internationale
  puis adminstratif, organisation.
*/





