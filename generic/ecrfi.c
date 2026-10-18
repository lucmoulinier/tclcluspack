#include "mainClusPack.h"

/**********************************************************************************/
/*                                                                                */
/*Procedure d'ecriture du fichier.dst des distances en vue du calcul de phylogenie*/
/*                                                                                */
/**********************************************************************************/

int ecritureFichierDistancesSequences(char *nomFichier,int nbIndividus,
				       individu_t **individus)
{
  /*declaration des variables*/
  FILE *fichier;
  int i,j;
  /*fin declaration des variables*/

  fichier=fopen(nomFichier,"w+");
  if(fichier==NULL)
    {
      return ECHEC;
    }

  fprintf(fichier,"%d\n",nbIndividus+1);
  for(i=0;i<nbIndividus;i++)
    {
      fprintf(fichier,"%d ",individus[i]->id);
      for(j=0;j<nbIndividus;j++)
	{
	  fprintf(fichier,"%f ",(float)individus[i]->valeursTraitees[individus[j]->id]);
	  if((j+1)%8==0)
	    {
	      fprintf(fichier,"\n");
	    }
	}
      fprintf(fichier,"1.0 \n");
    }
  fprintf(fichier,"fictive_sequence ");
  for(j=0;j<nbIndividus;j++)
    {
      fprintf(fichier,"1.0 ");
      if((j+1)%8==0)
	{
	  fprintf(fichier,"\n");
	}
    }
  fprintf(fichier,"0.0 \n");
  fclose(fichier);
  return SUCCES;
}

/******************************************************************/
/*                                                                */
/*Procedure d'ecriture du fichier resultat presentant les clusters*/
/*                                                                */
/******************************************************************/

void ecritureFichierClusters(int nbIndividus,int nbDimensions,individu_t **individus,
			       int nbClusters,int nbClustersSelection,
			       int writeCoordinates,char *nomFichier)
{
  /*declaration des variables*/
  int i,j,k,compteur,*taillesClusters,nbIndividusTotal,nbClustersRectifie;
  individu_t ***clusters;
  FILE *file;
  /*fin declaration des variables*/

  /*allocation memoire*/
  taillesClusters=(int *)ckalloc(sizeof(int)*(nbClusters+1));
  clusters=(individu_t ***)ckalloc(sizeof(individu_t *)*(nbClusters+1));
  /*fin allocation memoire*/

  nbIndividusTotal=nbIndividus; 
  printf("nbIndividusTotal %d\n",nbIndividusTotal);
 
  for(i=0;i<nbIndividus;i++) {
      nbIndividusTotal+=individus[i]->nbIndividusSimilaires;
    }

  printf("nbIndividusTotal %d\n",nbIndividusTotal);
  printf("nbClusters dans ecrfi : %d\n",nbClusters);
  for(i=0;i<nbClusters+1;i++)
    {
      clusters[i]=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividusTotal);
    }
  for(i=0;i<nbClusters+1;i++)
    {
      taillesClusters[i]=0;
    }

  for(i=0;i<nbIndividus;i++)
    {
      if(individus[i]->cluster>-1)
	{
	  clusters[individus[i]->cluster][taillesClusters[individus[i]->cluster]]=individus[i];
	  (taillesClusters[individus[i]->cluster])++;
	  
	  for(j=0;j<individus[i]->nbIndividusSimilaires;j++) {
	    clusters[individus[i]->cluster][taillesClusters[individus[i]->cluster]]=
	      individus[i]->individusSimilaires[j];
	    (taillesClusters[individus[i]->cluster])++;
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

  file=fopen(nomFichier,"w");
  fprintf(file,"Number of clusters : %d\n",nbClustersRectifie);
  
  compteur=0;
  for(i=0;i<nbClusters;i++)
    {
      if(taillesClusters[i]>0) {
	fprintf(file,"\nCluster %d ; size=%d\n",compteur,taillesClusters[i]);
	for(j=0;j<taillesClusters[i];j++)
	  {
	    fprintf(file,"%s",clusters[i][j]->nom);
	    if(writeCoordinates==OUI)
	      {
		for(k=0;k<nbDimensions;k++)
		  {
		    fprintf(file,"\t%f",clusters[i][j]->valeursBrutes[k]);
		  }
	      }	  
	    if(clusters[i][j]->description!=NULL)
	      {
		fprintf(file,"%s",clusters[i][j]->description);
	      }
	    else
	      {
		fprintf(file,"\n");
	      }
	  }
	compteur++;
      }
    }
  if(taillesClusters[nbClusters]>0)
    {
      fprintf(file,"\nunclustered points %d ; size=%d\n",i,taillesClusters[nbClusters]);
      for(j=0;j<taillesClusters[nbClusters];j++)
	{
	  fprintf(file,"%s",clusters[nbClusters][j]->nom);
	  if(writeCoordinates==OUI)
	    {
	      for(k=0;k<nbDimensions;k++)
		{
		  fprintf(file,"\t%f",clusters[nbClusters][j]->valeursTraitees[k]);
		}
	    }
	  if(clusters[i][j]->description!=NULL)
	    {
	      fprintf(file,"%s",clusters[i][j]->description);
	    }
	  else
	    {
	      fprintf(file,"\n");
	    }
	}
    }
  fclose(file);
 
  /*desallocation memoire*/
  ckfree(taillesClusters);
  for(i=0;i<nbClusters;i++)
    {
      ckfree(clusters[i]);
    }
  ckfree(clusters);
  /*fin desallocation memoire*/
}

/**********************************************************/
/*                                                        */
/*Procedure d'ecriture du fichier de l'alignement multiple*/
/*avec les sequences reordonnees suivant leurs groupes    */
/*                                                        */
/**********************************************************/

void ecritureFichierClustersAlignement(char *nomFichier,int longueurAlignement,
				       int nbIndividus,individu_t **individus,
				       char **sequences,int nbClusters)
{
  /*declaration des variables*/
  int i,j,k,*taillesClusters,nbIndividusTotal;
  individu_t ***clusters;
  FILE *file;
  /*fin declaration des variables*/

  /*allocation memoire*/
  taillesClusters=(int *)ckalloc(sizeof(int)*(nbClusters+1));
  clusters=(individu_t ***)ckalloc(sizeof(individu_t *)*(nbClusters+1));
  /*fin allocation memoire*/

  nbIndividusTotal=nbIndividus;
  for(i=0;i<nbClusters+1;i++)
    {
      clusters[i]=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividusTotal);
    }
  for(i=0;i<nbClusters+1;i++)
    {
      taillesClusters[i]=0;
    }

  for(i=0;i<nbIndividus;i++)
    {
      if(individus[i]->cluster>-1)
	{
	  clusters[individus[i]->cluster][taillesClusters[individus[i]->cluster]]=individus[i];
	  (taillesClusters[individus[i]->cluster])++;
	}
      else
	{
	  clusters[nbClusters][taillesClusters[nbClusters]]=individus[i];
	  (taillesClusters[nbClusters])++;
	}
    }

  file=fopen(nomFichier,"w");
  for(i=0;i<nbClusters;i++)
    {
      fprintf(file,">GROUP_%d\n",i+1);
      for(j=0;j<longueurAlignement;j++)
	{
	  fprintf(file,"-");
	}
      fprintf(file,"\n");
      for(j=0;j<taillesClusters[i];j++)
	{
	  fprintf(file,">%s\n",clusters[i][j]->nom);
	  for(k=0;k<longueurAlignement;k++)
	    {
	      if(sequences[clusters[i][j]->id][k]=='O')
		{
		  fprintf(file,"-");
		}
	      else
		{
		  fprintf(file,"%c",sequences[clusters[i][j]->id][k]);
		}
	    }
	  fprintf(file,"\n");
	}
    }
  if(taillesClusters[nbClusters]>0)
    {
      fprintf(file,">UNCLUSTERED\n");
      for(j=0;j<taillesClusters[nbClusters];j++)
	{
	  fprintf(file,">%s\n",clusters[nbClusters][j]->nom);
	  for(k=0;k<longueurAlignement;k++)
	    {
	      if(sequences[clusters[nbClusters][j]->id][k]=='O')
		{
		  fprintf(file,"-");
		}
	      else
		{
		  fprintf(file,"%c",sequences[clusters[nbClusters][j]->id][k]);
		}
	    }
	  fprintf(file,"\n");
	}
    }
  fclose(file);
 
  /*desallocation memoire*/
  ckfree(taillesClusters);
  for(i=0;i<nbClusters;i++)
    {
      ckfree(clusters[i]);
    }
  ckfree(clusters);
  /*fin desallocation memoire*/
}


/***************************************************/
/*                                                 */
/*Procedure d'ecriture des etiquettes d'un ensemble*/
/*                                                 */
/***************************************************/

void ecritureFichierEtiquettes(int tailleTestSet,individu_t **testSet,char *nomFichier) {
	/*declaration de variables*/
	int i;
	FILE *fichier;
	/*fin declaration de variables*/
	printf("fichier : %s\n",nomFichier);
	fichier=fopen(nomFichier,"w");
	for(i=0;i<tailleTestSet;i++) {
	  fprintf(fichier,"%s %d\n",testSet[i]->nom,testSet[i]->cluster);	
	}
	fclose(fichier);
}
