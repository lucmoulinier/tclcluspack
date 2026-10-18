#include "mainClusPack.h"
#include "tcl.h"
#include "TclDataOut.h"

/******************************************************************/
/*                                                                */
/*Procedure d'ecriture du fichier resultat presentant les clusters*/
/*                                                                */
/******************************************************************/

int RetourneClusters(Tcl_Interp *interp,
		     int nbIndividus,
		     int nbDimensions,
		     individu_t **individus,
		     int nbClusters,
		     int nbClustersSelection,
		     int writeCoordinates)
{
  /*declaration des variables*/
  int i,j,k,compteur,*taillesClusters,nbIndividusTotal,nbClustersRectifie, c;
  individu_t ***clusters;
  Tcl_Obj *ResultPtr, *IndPtr, *ClusPtr;
  /*fin declaration des variables*/

  /*allocation memoire*/
  taillesClusters = (int *)ckalloc(sizeof(int)*(nbClusters+1));
  clusters = (individu_t ***)ckalloc(sizeof(individu_t *)*(nbClusters+1));
  /*fin allocation memoire*/

  nbIndividusTotal = nbIndividus; 
  for(i=0;i<nbIndividus;i++) {
    nbIndividusTotal+=individus[i]->nbIndividusSimilaires;
  }
  
  for (i=0;i<nbClusters+1;i++) {
    clusters[i]=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividusTotal);
  }
  for(i=0;i<nbClusters+1;i++) {
    taillesClusters[i]=0;
  }
  
  //fprintf(stderr,"main loop\n");
  //fprintf(stderr,"nbIndividus %d\n",nbIndividus);
  //fflush(stderr);
  if(nbClusters == 0 && nbClusters != nbClustersSelection) {
    for(i=0;i<nbIndividus;i++)
      {
	individus[i]->cluster=0;
      }
  }

  for(i=0;i<nbIndividus;i++)
    {
      c = (int)individus[i]->cluster; 
      if(c > -1)
	{
	  //fprintf(stderr,"c %d taillesClusters %d\n",c,taillesClusters[c]);
	  //fflush(stderr);
	  clusters[c][taillesClusters[c]]=individus[i];
	  taillesClusters[c]++;
	  for(j=0;j<individus[i]->nbIndividusSimilaires;j++) {
	    clusters[c][taillesClusters[c]] = individus[i]->individusSimilaires[j];
	    taillesClusters[c]++;
	  }
	}
      else
	{
	  clusters[nbClusters][taillesClusters[nbClusters]]=individus[i];
	  (taillesClusters[nbClusters])++;
	  
	  for(j=0;j<individus[i]->nbIndividusSimilaires;j++) {
	    clusters[nbClusters][taillesClusters[nbClusters]]=
	      individus[i]->individusSimilaires[j];
	    (taillesClusters[nbClusters])++;
	  }
	  
	}
    }

  /*calcul du vrai nombre de clusters, en supprimant les clusters vides*/
  nbClustersRectifie=0;
  for(i=0;i<nbClusters+1;i++) {
    if(taillesClusters[i]>0) {
      nbClustersRectifie++;
    }
  }
  ResultPtr = Tcl_NewListObj(0,NULL);
  compteur=0;
  for(i=0;i<nbClusters;i++)
    {
      if (taillesClusters[i]>0) {
	ClusPtr = Tcl_NewListObj(0,NULL);
	for (j=0;j<taillesClusters[i];j++)
	  {
	    IndPtr = Tcl_NewListObj(0,NULL);
	    Tcl_ListObjAppendElement(interp,IndPtr,Tcl_NewStringObj(clusters[i][j]->nom,-1));
	    if (writeCoordinates==OUI)
	      {
		for (k=0;k<nbDimensions;k++)
		  {
		    Tcl_ListObjAppendElement(interp,IndPtr,Tcl_NewDoubleObj(clusters[i][j]->valeursBrutes[k]));
		  }
	      }	  
	    if(clusters[i][j]->description!=NULL)
	      {
		//fprintf(file,"%s",clusters[i][j]->description);
	      }
	    else
	      {
		//fprintf(file,"\n");
	      }
	    Tcl_ListObjAppendElement(interp,ClusPtr,IndPtr);
	  }
	compteur++;
	Tcl_ListObjAppendElement(interp,ResultPtr,ClusPtr);
      }
    }
  if (taillesClusters[nbClusters]>0) {
    ClusPtr = Tcl_NewListObj(0,NULL);
    for (j=0;j<taillesClusters[nbClusters];j++) {
      IndPtr = Tcl_NewListObj(0,NULL);
      Tcl_ListObjAppendElement(interp,IndPtr,Tcl_NewStringObj(clusters[nbClusters][j]->nom,-1));
      if (writeCoordinates==OUI) {
	for(k=0;k<nbDimensions;k++) {
	  Tcl_ListObjAppendElement(interp,IndPtr,Tcl_NewDoubleObj(clusters[nbClusters][j]->valeursBrutes[k]));
	}
      }
      Tcl_ListObjAppendElement(interp,ClusPtr,IndPtr);
    }
    Tcl_ListObjAppendElement(interp,ResultPtr,ClusPtr);
  }
  //fclose(file);
 
  /*desallocation memoire*/
  ckfree((char *)taillesClusters);
  for(i=0;i<nbClusters+1;i++) {
    ckfree((char *)clusters[i]);
  }
  ckfree((char *)clusters);
  /*fin desallocation memoire*/

  Tcl_SetObjResult(interp,ResultPtr);
  return TCL_OK;
}


