#include "mainClusPack.h"
#include "matrice.h"
#include "tcl.h"

/********************************************/
/*                                          */
/*Procedure de calcul de la densite uniforme*/
/*                                          */
/********************************************/

double calculLogDensiteUniformeGaussienne(int nbIndividus,int nbDimensions,individu_t **individus) {
  /*declaration de variables*/
  int i,j;
  double *valeursMinimum,*valeursMaximum,logDensiteLoiUniforme;
  /*fin declaration de variables*/

  /*desallocation memoire*/
  valeursMinimum=(double *)ckalloc(sizeof(double)*nbDimensions);
  valeursMaximum=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/

  for(i=0;i<nbDimensions;i++)
    {
      valeursMinimum[i]=individus[0]->valeursTraitees[i];
      valeursMaximum[i]=individus[0]->valeursTraitees[i];
    }
  
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  if(valeursMinimum[j]>individus[i]->valeursTraitees[j])
	    {
	      valeursMinimum[j]=individus[i]->valeursTraitees[j];
	    }
	  if(valeursMaximum[j]<individus[i]->valeursTraitees[j])
	    {
	      valeursMaximum[j]=individus[i]->valeursTraitees[j];
	    }
	}
    }
  logDensiteLoiUniforme=0.0;
  for(i=0;i<nbDimensions;i++)
    {
      /*une dimension constante donnerait log(0)*/
      if(valeursMaximum[i]>valeursMinimum[i]) {
	logDensiteLoiUniforme-=log(3.0*(valeursMaximum[i]-valeursMinimum[i]));
      }
    }

  /*desallocation memoire*/
  ckfree(valeursMinimum);
  ckfree(valeursMaximum);
  /*fin desallocation memoire*/

  return logDensiteLoiUniforme;
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calculVraisemblancesPointsGaussienne(int nbIndividus,
					  int nbDimensions,
					  individu_t **individus,
					  int nbClusters,
					  cluster_t *clusters,
					  int diag)
{
  /*declaration de variables*/
  int i,j,k;
  double **sigma,*sigmaDiag,norme;
  double logDeterminant,*difference,commun;
  /*fin declaration de variables*/

  //printf("dimension : %d\n",nbDimensions);

  /*allocation memoire*/
  if(diag==NON) {
    sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
    for(i=0;i<nbDimensions;i++)
      {
	sigma[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
      }
  } else {
    sigmaDiag=(double *)ckalloc(sizeof(double)*nbDimensions);
  }  
  difference=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/
  
  for(i=0;i<nbClusters;i++)
    {
      /*on inverse la matrice de covariance et on recupere au passage*/
      /*le log du determinant de cette meme matrice (de covariance)*/
      if(diag==NON) {
	/*on inverse une copie de la matrice de covariance dans sigma
	  (regularisee si elle est singuliere)*/
	logDeterminant=inversionMatriceRegularisee(nbDimensions,clusters[i].sigma,sigma);
      } else {
	/*on copie la matrice de covariance dans sigma*/
	for(j=0;j<nbDimensions;j++)
	  {
	    sigmaDiag[j]=clusters[i].sigmaDiag[j];
	    //printf("sigma : %f\n",sigmaDiag[j]);
	    if(sigmaDiag[j]==0) {
	      sigmaDiag[j]=1.0;
	    }
	  }
	logDeterminant=0;
	for(j=0;j<nbDimensions;j++) {
	  logDeterminant+=log(sigmaDiag[j]);
	  sigmaDiag[j]=1.0/sigmaDiag[j];
	}
      }
      //  printf("log_det : %f\n",logDeterminant);
      commun=-logDeterminant*0.5-(double)nbDimensions*0.5*log(2.0*M_PI)+
	log(clusters[i].tau);
      //log((double)clusters[i].nbIndividus/(double)nbIndividus);
      
      /*calcul des log-vraisemblances*/
      for(j=0;j<nbIndividus;j++)
	{
	  for(k=0;k<nbDimensions;k++)
	    {
	      difference[k]=individus[j]->valeursTraitees[k]-clusters[i].mu[k];
	    }
	  
	  /*individus[j]->vrais[i]=-log_determinant/2.0-(double)nbDimensions/2.0*log(2.0*M_PI);*/
	  if(diag==NON) {
	    norme=calculNorme(nbDimensions,difference,sigma);
	  } else {
	    norme=0;
	    for(k=0;k<nbDimensions;k++) {
	      norme+=difference[k]*difference[k]*sigmaDiag[k];
	    }
	  }
	  individus[j]->vrais[i]=commun-norme*0.5;
	  /*individus[j]->vrais[i]+=log((double)clusters[i].nbIndividus/(double)nbIndividus);*/
	  //printf("vrais indi %f\n",individus[j]->vrais[i]);
	  
	}
    }
  
  /*desallocation memoire*/
  if(diag==NON) {
    for(i=0;i<nbDimensions;i++)
      {
	ckfree((char *)sigma[i]);
      }
    ckfree((char *)sigma);
  } else {
    ckfree((char *)sigmaDiag);
  }
  ckfree((char *)difference);
  /*fin desallocation memoire*/
}

/**************************************************/
/*                                                */
/*Procedure de calcul du mu moyen pour initialiser*/
/*la maximisation de la vraisemblance             */
/*                                                */
/**************************************************/

void calculMuMoyen(int nbPoints,double **coordonnees,int nbDimensions,
		     double *mu)
{
  /*declaration de variables*/
  int i,j;
  /*fin declaration de variables*/

  for(j=0;j<nbDimensions;j++)
    {
      mu[j]=0;
    }

  for(i=0;i<nbPoints;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  mu[j]+=coordonnees[i][j];
	}
    }
  for(j=0;j<nbDimensions;j++)
    {
      mu[j]/=(double)nbPoints;
    }
}


/****************************************************************************/
/*                                                                          */
/*Procedure de maximisation de la vraisemblance pour des donnees gaussiennes*/
/*                                                                          */
/****************************************************************************/
  
int maximisationVraisemblanceGaussienne(double **points,int nbDimensions,cluster_t cluster,
					int diag) {
  
  /*calcul du mu moyen*/
  calculMuMoyen(cluster.nbIndividus,points,nbDimensions,cluster.mu);
  
  /*estimation de la matrice de covariance*/
  if(diag==OUI) {
    calculMatriceVariance(nbDimensions,cluster.nbIndividus,points,cluster.sigmaDiag);
  } else {
    calculMatriceCovariance(nbDimensions,cluster.nbIndividus,points,cluster.sigma);
  }
  /*  int i;
  for(i=0;i<nbDimensions;i++) {
    printf("sigmaDiag i %f\n",cluster.sigmaDiag[i]);
    }*/
  return OK;
}
