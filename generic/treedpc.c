#include "mainClusPack.h"
#include "dpc.h"
#include "treecut.h"
#include "tcl.h"

/********************************************************************************/
/*                                                                              */
/*Procedure de decoupage d'arbre avec le test anova sur les densites de presence*/
/*                                                                              */
/********************************************************************************/

void decoupageArbreDPC(int nbIndividus,
		       int nbDimensions,
		       individu_t **individus,
		       noeud_t *racine,
		       int *nbClusters,
		       int typeDonnees,
		       int nbClustersSelection,
		       int donneesNormalisees,
		       int typeDensite,
		       int nbVoisins)
{
  /*declaration de variables*/
  int i,j,auMoinsUnNoeudSecable,nbNoeudsActifs,test,nbFoisQuonSepare;
  int nbSimulations,noeudCourant,compteur;
  double **coordonneesIndividus,*cotesHyperpave,*valeursMinimum,*valeursMaximum,risque;
  double qualite;
  noeud_t **noeuds,**feuilles;
  cluster_t groupePere,groupe1,groupe2;
  /*fin declaration de variables*/

  /*allocation memoire*/
  noeuds=(noeud_t **)ckalloc(sizeof(noeud_t *)*nbIndividus);
  feuilles=(noeud_t **)ckalloc(sizeof(noeud_t *)*nbIndividus);
  groupePere.individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  groupe1.individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  groupe2.individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
  coordonneesIndividus=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      coordonneesIndividus[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  cotesHyperpave=(double *)ckalloc(sizeof(double)*nbDimensions);
  valeursMinimum=(double *)ckalloc(sizeof(double)*nbDimensions);
  valeursMaximum=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  /*on recupere toutes les feuilles*/
  compteur=0;
  chercheFeuilles(racine,&compteur,feuilles);
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  coordonneesIndividus[i][j]=feuilles[i]->individu->valeursTraitees[j];
	}
    }

  qualite=0;
  risque=0.01;
  if(typeDonnees==ALIGNEMENT)
    {
      risque=0.01;
    }
  else if(typeDensite==DENSITE2)
    {
      risque=0.01;
    }
  
  nbSimulations=10;
  if(typeDonnees==ALIGNEMENT)
    { 
      for(i=0;i<nbDimensions;i++)
	{
	  valeursMinimum[i]=0;
	  cotesHyperpave[i]=1.0;
	}
    }
  else
    {
      /*on cherche les valeurs minimum et maximum pour chaque axe*/
      for(i=0;i<nbDimensions;i++)
	{
	  valeursMinimum[i]=coordonneesIndividus[0][i];
	  valeursMaximum[i]=coordonneesIndividus[0][i];
	}
      
      for(i=0;i<nbIndividus;i++)
	{
	  for(j=0;j<nbDimensions;j++)
	    {
	      if(valeursMinimum[j]>coordonneesIndividus[i][j])
		{
		  valeursMinimum[j]=coordonneesIndividus[i][j];
		}
	      if(valeursMaximum[j]<coordonneesIndividus[i][j])
		{
		  valeursMaximum[j]=coordonneesIndividus[i][j];
		}
	    }
	}
      
      /*on calcule la longueur des cotes de l'hyperpave*/
      for(i=0;i<nbDimensions;i++)
	{
	  cotesHyperpave[i]=valeursMaximum[i]-valeursMinimum[i];
	}
    }
    
  nbNoeudsActifs=1;
  noeuds[0]=racine;
  noeuds[0]->secable=OUI;
  *nbClusters=0;
  while(1)
    {
      auMoinsUnNoeudSecable=NON;
      for(i=0;i<nbNoeudsActifs;i++)
	{
	  if(noeuds[i]->secable==OUI)
	    {
	      auMoinsUnNoeudSecable=OUI;
	      noeudCourant=i;
	      break;
	    }
	}
      if(auMoinsUnNoeudSecable==NON)
	{
	  break;
	}
      else
	{ 
	  /*on cherche les feuilles de la branche du noeud i*/
	  groupePere.nbIndividus=0;
	  chercheFeuilles(noeuds[noeudCourant],&(groupePere.nbIndividus),feuilles);
	  
	  for(i=0;i<groupePere.nbIndividus;i++)
	    {
	      groupePere.individus[i]=feuilles[i]->individu;
	    }

	  if(groupePere.nbIndividus<=2)
	    {
	      noeuds[noeudCourant]->secable=NON;
	      for(j=0;j<groupePere.nbIndividus;j++)
		{
		  groupePere.individus[j]->cluster=*nbClusters;
		}
	      (*nbClusters)++;
	    }
	  else
	    {
	      /*on cherche les feuilles de la branche 1 du noeud i*/
	      groupe1.nbIndividus=0;
	      chercheFeuilles(noeuds[noeudCourant]->copain1,&(groupe1.nbIndividus),feuilles);
	      
	      for(i=0;i<groupe1.nbIndividus;i++)
		{
		  groupe1.individus[i]=feuilles[i]->individu;
		}
	      
	      /*on cherche les feuilles de la branche 2 du noeud i*/
	      groupe2.nbIndividus=0;
	      chercheFeuilles(noeuds[noeudCourant]->copain2,&(groupe2.nbIndividus),feuilles);
	      
	      for(i=0;i<groupe2.nbIndividus;i++)
		{
		  groupe2.individus[i]=feuilles[i]->individu;
		}
	      
	      if(nbClustersSelection==DPC)
		{
		  /*test de division possible par DPC*/
		  nbFoisQuonSepare=0;
		  for(j=0;j<nbSimulations;j++)
		    {
		      /*on test si la separation de groupePere en groupe1 et
			groupe2 est pertinente*/
		      test=testDPC(nbIndividus,nbDimensions,&groupePere,&groupe1,&groupe2,
				   coordonneesIndividus,cotesHyperpave,valeursMinimum,risque,
				   typeDonnees,typeDensite,donneesNormalisees,nbSimulations);
		      if(test==OUI)
			{
			  nbFoisQuonSepare++;
			}
		    }
		  if(nbFoisQuonSepare>nbSimulations/2)
		    {
		      test=OUI;
		      qualite+=(double)nbFoisQuonSepare/(double)nbSimulations;
		    }
		  else
		    {
		      test=NON;
		      qualite+=(double)(nbSimulations-nbFoisQuonSepare)/(double)nbSimulations;
		    }
		}
	      
	      if(test==OUI)
		{
		  noeuds[nbNoeudsActifs]=noeuds[noeudCourant]->copain1;
		  noeuds[nbNoeudsActifs+1]=noeuds[noeudCourant]->copain2;
		  nbNoeudsActifs+=2;
		  noeuds[noeudCourant]->secable=NON;
		  noeuds[noeudCourant]->copain1->secable=OUI;
		  noeuds[noeudCourant]->copain2->secable=OUI;
		}
	      else
		{
		  for(j=0;j<groupePere.nbIndividus;j++)
		    {
		      groupePere.individus[j]->cluster=*nbClusters;
		    }
		  (*nbClusters)++;
		  noeuds[noeudCourant]->secable=NON;
		}
	    }
	}
    }

  /*desallocation memoire*/
  ckfree((char *)noeuds);
  ckfree((char *)feuilles);
  ckfree((char *)groupePere.individus);
  ckfree((char *)groupe1.individus);
  ckfree((char *)groupe2.individus); 
  for(i=0;i<nbIndividus;i++)
    {
      ckfree((char *)coordonneesIndividus[i]);
    }
  ckfree((char *)coordonneesIndividus);
  ckfree((char *)valeursMinimum); 
  ckfree((char *)valeursMaximum);
  ckfree((char *)cotesHyperpave);
  /*fin desallocation memoire*/
}
