#include "mainClusPack.h"
#include "kernelKmeans.h"
#include "tools.h"
#include "kernels.h"

/****************************************************/
/*                                                  */
/*Procedure de realisation de clustering par k-means*/
/*en partant d'un debut de clustering               */ 
/*                                                  */
/****************************************************/
/*typeNoyau peut valoir AUCUN_NOYAU ou RADIAL_GAUSSIAN*/
/*2026-10-08 : corrections
  - germes tires par floor(n-1)*runif : le dernier individu n'etait jamais
    tire et la boucle de tirage ne terminait pas si nbClusters >= nbIndividus.
    -> initialisation k-means++ dans l'espace du noyau (germe suivant tire
    avec une probabilite proportionnelle a sa distance au carre au germe le
    plus proche), qui evite de placer deux germes dans le meme nuage ;
  - au debut de chaque relance, la meilleure solution etait remplacee sans
    comparaison par l'affectation initiale de cette relance, et la solution
    finale d'une relance n'etait pas enregistree quand elle avait converge
    -> la meilleure inertie est comparee apres chaque affectation ;
  - arret quand plus aucun individu ne change de groupe (au lieu d'un ecart
    d'inertie), 10 relances au lieu de 2, 100 iterations au maximum.*/

/*distance au carre entre deux individus dans l'espace du noyau*/
static double distanceNoyau(int nbDimensions,individu_t *point1,individu_t *point2,
			    int choixNoyau,int accelerationCalculs,double **kernels) {
  if(accelerationCalculs==NON) {
    return kernel(nbDimensions,point1,point1,choixNoyau)+
      kernel(nbDimensions,point2,point2,choixNoyau)-
      2.0*kernel(nbDimensions,point1,point2,choixNoyau);
  }
  return kernels[point1->id][point1->id]+kernels[point2->id][point2->id]-
    2.0*kernels[point1->id][point2->id];
}

/*initialisation k-means++ : germes dans pointsSelectionnes[0..nbClusters-1]*/
static void choixGermesKernelKmeans(int nbIndividus,int nbDimensions,individu_t **individus,
				    int nbClusters,int *pointsSelectionnes,double *distancesMin,
				    int choixNoyau,int accelerationCalculs,double **kernels) {
  int i,j;
  double somme,seuil,d;

  pointsSelectionnes[0]=(int)floor((double)nbIndividus*runif(0,1));
  if(pointsSelectionnes[0]>=nbIndividus) {
    pointsSelectionnes[0]=nbIndividus-1;
  }
  for(i=0;i<nbIndividus;i++) {
    distancesMin[i]=distanceNoyau(nbDimensions,individus[i],individus[pointsSelectionnes[0]],
				  choixNoyau,accelerationCalculs,kernels);
  }
  for(j=1;j<nbClusters;j++) {
    somme=0;
    for(i=0;i<nbIndividus;i++) {
      if(distancesMin[i]>0) {
	somme+=distancesMin[i];
      }
    }
    pointsSelectionnes[j]=-1;
    if(somme>0) {
      seuil=somme*runif(0,1);
      for(i=0;i<nbIndividus;i++) {
	if(distancesMin[i]>0) {
	  pointsSelectionnes[j]=i;
	  seuil-=distancesMin[i];
	  if(seuil<0) {
	    break;
	  }
	}
      }
    } else {
      /*individus restants confondus avec les germes : premier non choisi*/
      for(i=0;i<nbIndividus;i++) {
	int k;
	for(k=0;k<j;k++) {
	  if(pointsSelectionnes[k]==i) break;
	}
	if(k==j) {
	  pointsSelectionnes[j]=i;
	  break;
	}
      }
    }
    for(i=0;i<nbIndividus;i++) {
      d=distanceNoyau(nbDimensions,individus[i],individus[pointsSelectionnes[j]],
		      choixNoyau,accelerationCalculs,kernels);
      if(d<distancesMin[i]) {
	distancesMin[i]=d;
      }
    }
    distancesMin[pointsSelectionnes[j]]=0;
  }
}

int kernelKmeans(int nbIndividus,
		 int nbDimensions,
		 individu_t **individus,
		 int nbClusters,
		 int choixNoyau,
		 int accelerationCalculs,
		 double **kernels) {
  /*declaration des variables*/
  int i,*solution,*affectationPrecedente,compteurIterations1,compteurIterations2;
  int nbKmeans,nbIterationsMax,*pointsSelectionnes,resultat,changement;
  double inertieIntraClasses,meilleureInertieIntraClasses,*distancesMin;
  cluster_t *clusters;
  /*fin eclaration des variables*/

  if((nbClusters<1)||(nbClusters>nbIndividus)) {
    return ECHEC;
  }

  /*allocation memoire*/
  clusters=(cluster_t *)ckalloc(sizeof(cluster_t)*nbClusters);
  for(i=0;i<nbClusters;i++)
    {
      clusters[i].individus=(individu_t **)ckalloc(sizeof(individu_t *)*nbIndividus);
    }
  solution=(int *)ckalloc(sizeof(int)*nbIndividus);
  affectationPrecedente=(int *)ckalloc(sizeof(int)*nbIndividus);
  pointsSelectionnes=(int *)ckalloc(sizeof(int)*nbClusters);
  distancesMin=(double *)ckalloc(sizeof(double)*nbIndividus);
  /*fin allocation memoire*/

  nbKmeans=10;
  nbIterationsMax=100;
  meilleureInertieIntraClasses=DBL_MAX;
  resultat=ECHEC;

  /*kmeans*/
  for(compteurIterations1=0;compteurIterations1<nbKmeans;compteurIterations1++)
    {
      choixGermesKernelKmeans(nbIndividus,nbDimensions,individus,nbClusters,
			      pointsSelectionnes,distancesMin,
			      choixNoyau,accelerationCalculs,kernels);
      for(i=0;i<nbClusters;i++)
	{
	  clusters[i].nbIndividus=1;
	  clusters[i].individus[0]=individus[pointsSelectionnes[i]];
	}

      /*on affecte les individus a leurs groupes les plus proches*/
      if(affectationIndividusKernelKmeans(nbDimensions,nbIndividus,nbClusters,clusters,
					  individus,&inertieIntraClasses,1,
					  choixNoyau,accelerationCalculs,kernels)==ECHEC)
	{
	  continue;
	}
      for(compteurIterations2=0;compteurIterations2<nbIterationsMax;compteurIterations2++)
	{
	  for(i=0;i<nbIndividus;i++)
	    {
	      affectationPrecedente[i]=individus[i]->cluster;
	    }
	  /*on affecte les individus a leurs groupes les plus proches*/
	  if(affectationIndividusKernelKmeans(nbDimensions,nbIndividus,nbClusters,clusters,
					      individus,&inertieIntraClasses,0,
					      choixNoyau,accelerationCalculs,kernels)==ECHEC)
	    {
	      break;
	    }
	  changement=0;
	  for(i=0;i<nbIndividus;i++)
	    {
	      if(affectationPrecedente[i]!=individus[i]->cluster)
		{
		  changement=1;
		  break;
		}
	    }
	  if(!changement)
	    {
	      break;
	    }
	}
      if(inertieIntraClasses<meilleureInertieIntraClasses)
	{
	  for(i=0;i<nbIndividus;i++)
	    {
	      solution[i]=individus[i]->cluster;
	    }
	  meilleureInertieIntraClasses=inertieIntraClasses;
	  resultat=SUCCES;
	}
    }
  if(resultat==SUCCES)
    {
      for(i=0;i<nbIndividus;i++)
	{
	  individus[i]->cluster=solution[i];
	}
    }
  
  /*desallocation memoire*/
  for(i=0;i<nbClusters;i++)
    {
      ckfree(clusters[i].individus);
    }
  ckfree(clusters);
  ckfree(solution);
  ckfree(affectationPrecedente);
  ckfree(pointsSelectionnes);
  ckfree(distancesMin);
  /*fin desallocation memoire*/

  return resultat;
}

/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*pour la methode de kernel kmeans                 */
/*                                                 */
/***************************************************/

int affectationIndividusKernelKmeans(int nbDimensions,
				     int nbIndividus,
				     int nbClusters,
				     cluster_t *clusters,
				     individu_t **individus,
				     double *inertieIntraClasses,
				     int initialisation,
				     int choixNoyau,
				     int accelerationCalculs,
				     double **kernels)
{
  /*declaration des variables*/
  int i,j,k,clusterLePlusProche,*nbIndividusParClasse,**classes,resultat=SUCCES;
  double distanceCourante,plusPetiteDistance;
  /*fin declaration des variables*/

  /*allocation memoire*/
  nbIndividusParClasse=(int *)ckalloc(sizeof(int)*nbClusters);
  classes=(int **)ckalloc(sizeof(int *)*nbClusters);
  for(i=0;i<nbClusters;i++) {
    classes[i]=(int *)ckalloc(sizeof(int)*nbIndividus);
  }
  /*fin allocation memoire*/

  for(i=0;i<nbClusters;i++) {
      nbIndividusParClasse[i]=0;
    }

  if(initialisation==1) {
    if(accelerationCalculs==NON) {
      for(i=0;i<nbIndividus;i++) {
	plusPetiteDistance=DBL_MAX;
	clusterLePlusProche=-1;
	for(j=0;j<nbClusters;j++) {
	  distanceCourante=kernel(nbDimensions,individus[i],individus[i],choixNoyau)+
	    kernel(nbDimensions,clusters[j].individus[0],clusters[j].individus[0],choixNoyau)-
	    2.0*kernel(nbDimensions,individus[i],clusters[j].individus[0],choixNoyau);
	  if(distanceCourante<plusPetiteDistance) {
	    plusPetiteDistance=distanceCourante;
	    clusterLePlusProche=j;
	  }
	}
	if(clusterLePlusProche<0) {
	  /*distances non definies (NaN)*/
	  resultat=ECHEC;
	  break;
	}
	individus[i]->cluster=clusterLePlusProche;
	classes[clusterLePlusProche][nbIndividusParClasse[clusterLePlusProche]]=i;
	(nbIndividusParClasse[clusterLePlusProche])++;
      }
    } else {
      /*acceleration des calculs*/
      for(i=0;i<nbIndividus;i++) {
	plusPetiteDistance=DBL_MAX;
	clusterLePlusProche=-1;
	for(j=0;j<nbClusters;j++) {
	  distanceCourante=kernels[individus[i]->id][individus[i]->id]+
	    kernels[clusters[j].individus[0]->id][clusters[j].individus[0]->id]-
	    2.0*kernels[individus[i]->id][clusters[j].individus[0]->id];
	  if(distanceCourante<plusPetiteDistance) {
	    plusPetiteDistance=distanceCourante;
	    clusterLePlusProche=j;
	  }
	}
	if(clusterLePlusProche<0) {
	  /*distances non definies (NaN)*/
	  resultat=ECHEC;
	  break;
	}
	individus[i]->cluster=clusterLePlusProche;
	classes[clusterLePlusProche][nbIndividusParClasse[clusterLePlusProche]]=i;
	(nbIndividusParClasse[clusterLePlusProche])++;
      }
    }
  } else {
    if(accelerationCalculs==NON) {
      for(i=0;i<nbIndividus;i++)
	{
	  plusPetiteDistance=DBL_MAX;
	  clusterLePlusProche=-1;
	  for(j=0;j<nbClusters;j++)
	    {
	      if(clusters[j].nbIndividus==0) {
		/*groupe vide : division par zero (NaN) avant 2026-10-08*/
		continue;
	      }
	      distanceCourante=0;
	      for(k=0;k<clusters[j].nbIndividus;k++) {
		distanceCourante-=kernel(nbDimensions,individus[i],clusters[j].individus[k],choixNoyau);
	      }
	      distanceCourante*=2.0/(double)clusters[j].nbIndividus;
	      distanceCourante+=kernel(nbDimensions,individus[i],individus[i],choixNoyau)+
		clusters[j].sommeKxixj/(double)(clusters[j].nbIndividus*clusters[j].nbIndividus);
	      if(distanceCourante<plusPetiteDistance)
		{
		  plusPetiteDistance=distanceCourante;
		  clusterLePlusProche=j;
		}
	    }
	  if(clusterLePlusProche<0) {
	    /*distances non definies (NaN)*/
	    resultat=ECHEC;
	    break;
	  }
	  individus[i]->cluster=clusterLePlusProche;
	  classes[clusterLePlusProche][nbIndividusParClasse[clusterLePlusProche]]=i;
	  (nbIndividusParClasse[clusterLePlusProche])++;
	}
    } else {
      /*acceleration des calculs*/
      for(i=0;i<nbIndividus;i++) {
	  plusPetiteDistance=DBL_MAX;
	  clusterLePlusProche=-1;
	  for(j=0;j<nbClusters;j++)
	    {
	      if(clusters[j].nbIndividus==0) {
		/*groupe vide : division par zero (NaN) avant 2026-10-08*/
		continue;
	      }
	      distanceCourante=0;
	      for(k=0;k<clusters[j].nbIndividus;k++) {
		distanceCourante-=kernels[individus[i]->id][clusters[j].individus[k]->id];
	      }
	      distanceCourante*=2.0/(double)clusters[j].nbIndividus;
	      distanceCourante+=kernels[individus[i]->id][individus[i]->id]+
		clusters[j].sommeKxixj/(double)(clusters[j].nbIndividus*clusters[j].nbIndividus);
	      if(distanceCourante<plusPetiteDistance)
		{
		  plusPetiteDistance=distanceCourante;
		  clusterLePlusProche=j;
		}
	    }
	  if(clusterLePlusProche<0) {
	    /*distances non definies (NaN)*/
	    resultat=ECHEC;
	    break;
	  }
	  individus[i]->cluster=clusterLePlusProche;
	  classes[clusterLePlusProche][nbIndividusParClasse[clusterLePlusProche]]=i;
	  (nbIndividusParClasse[clusterLePlusProche])++;
	}
    }
  }


  if(resultat==ECHEC) {
    *inertieIntraClasses=DBL_MAX;
    goto liberation;
  }

  /*mise a jour des clusters*/
  for(i=0;i<nbClusters;i++) {
    clusters[i].nbIndividus=nbIndividusParClasse[i];
    for(j=0;j<clusters[i].nbIndividus;j++) {
      clusters[i].individus[j]=individus[classes[i][j]];
    }
  }

  *inertieIntraClasses=0;
  if(accelerationCalculs==NON) {
    for(i=0;i<nbClusters;i++)
      {
	clusters[i].sommeKxixj=0;
	for(j=0;j<clusters[i].nbIndividus;j++) {
	  *inertieIntraClasses+=kernel(nbDimensions,clusters[i].individus[j],clusters[i].individus[j],choixNoyau);
	  clusters[i].sommeKxixj+=kernel(nbDimensions,clusters[i].individus[j],clusters[i].individus[j],choixNoyau);
	  for(k=j+1;k<clusters[i].nbIndividus;k++) {
	    clusters[i].sommeKxixj+=
	      2.0*kernel(nbDimensions,clusters[i].individus[j],clusters[i].individus[k],choixNoyau);
	  }
	}
	if(clusters[i].nbIndividus>0) {
	  *inertieIntraClasses-=clusters[i].sommeKxixj/(double)clusters[i].nbIndividus;
	}
      }  
  } else {
    /*acceleration des calculs*/
    for(i=0;i<nbClusters;i++)
      {
	clusters[i].sommeKxixj=0;
	for(j=0;j<clusters[i].nbIndividus;j++) {
	  *inertieIntraClasses+=kernels[clusters[i].individus[j]->id][clusters[i].individus[j]->id];
	  clusters[i].sommeKxixj+=kernels[clusters[i].individus[j]->id][clusters[i].individus[j]->id];
	  for(k=j+1;k<clusters[i].nbIndividus;k++) {
	    clusters[i].sommeKxixj+=
	      2.0*kernels[clusters[i].individus[j]->id][clusters[i].individus[k]->id];
	  }
	}
	if(clusters[i].nbIndividus>0) {
	  *inertieIntraClasses-=clusters[i].sommeKxixj/(double)clusters[i].nbIndividus;
	}
      } 
  }

 liberation:
  /*desallocation memoire*/
  ckfree(nbIndividusParClasse);
  for(i=0;i<nbClusters;i++) {
    ckfree(classes[i]);
  }
  ckfree(classes);
  /*fin desallocation memoire*/

  return resultat;
}
