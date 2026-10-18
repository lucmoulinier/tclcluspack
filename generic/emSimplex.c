#include "mainClusPack.h"
#include "space.h"
#include "emSimplex.h"

/********************************************/
/*                                          */
/*Procedure de calcul de la densite uniforme*/
/*                                          */
/********************************************/

double calculLogDensiteUniformeSimplex(int nbIndividus,int nbDimensions,individu_t **individus) {
  /*declaration de variables*/
  double logDensiteLoiUniforme,volumeSimplex;
  /*fin declaration de variables*/
  
  /*calcul du volume d'un tetraedre*/
  calculVolumeSimplex(nbDimensions-1,&volumeSimplex);

  logDensiteLoiUniforme=-log(volumeSimplex);

  return logDensiteLoiUniforme;
}

/***********************************************/
/*                                             */
/*Procedure de maximisation de la vraisemblance*/
/*                                             */
/***********************************************/

int maximisationVraisemblanceSimplex(double **points,int nbDimensions,cluster_t cluster)
{
  /*declaration de variables*/
  int i,j,k,nbIterations,nbIterationsMax,probleme;
  double somme,sommeAlphas,*variances,*gradient,normeDeriveeCarre;
  double *nouveauxAlphas,logVraisemblanceCourante,nouvelleLogVraisemblance;
  double **inverseHessien,denominateur,a,alpha0Max,*fj;
  double differenceAlphas,sommeInverses;

  /*pour tests*/
  double premiersAlphas[50];
  /*fin declaration de variables*/

  probleme=NON;

  /*allocation memoire*/
  variances=(double *)ckalloc(sizeof(double)*nbDimensions);
  gradient=(double *)ckalloc(sizeof(double)*nbDimensions);
  inverseHessien=(double **)ckalloc(sizeof(double *)*nbDimensions);
  for(i=0;i<nbDimensions;i++)
    {
      inverseHessien[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  nouveauxAlphas=(double *)ckalloc(sizeof(double)*nbDimensions);
  fj=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*  cluster.alphas=(double *)ckalloc(sizeof(double)*nbDimensions);*/
  
  /*fin allocation memoire*/
  gradient[0]=5.0;
  nbIterationsMax=100;
	  
    /*on calcule les fj, estimation par les moments*/
  sommeInverses=0;
  for(k=0;k<nbDimensions;k++)
    {
      fj[k]=0;
      for(j=0;j<cluster.nbIndividus;j++)
	{
	  fj[k]+=points[j][k];
	}
      fj[k]/=(double)cluster.nbIndividus;
      sommeInverses+=1.0/fj[k];
    }    	  

  /*on calcule a en resolvant une equation du second degre apres approximation de la logvraisemblance*/
  a=0;
  for(k=0;k<nbDimensions;k++)
    {
      somme=0.0;
      for(j=0;j<cluster.nbIndividus;j++)
	{
	  somme+=log(points[j][k]);
	}
      a+=fj[k]*(somme-(double)cluster.nbIndividus*log(fj[k]));
    }
  
  alpha0Max=-(double)(cluster.nbIndividus*(nbDimensions-1))*0.5-
    sqrt((double)cluster.nbIndividus*((double)cluster.nbIndividus/4.0*
				       (double)(nbDimensions-1)*(double)(nbDimensions-1)-
				       a*(sommeInverses-1.0)/3.0));
  alpha0Max/=2.0*a;
  sommeAlphas=0;	  
  for(j=0;j<nbDimensions;j++)
    {
      cluster.alphas[j]=fj[j]*alpha0Max;
      sommeAlphas+=cluster.alphas[j];
      premiersAlphas[j]=cluster.alphas[j];
      if(!((premiersAlphas[j]>ZERO_LIMIT)&&(premiersAlphas[j]<999999999.0)))
	{
	  printf("%f ; %f ; %f ; %f\n",premiersAlphas[j],alpha0Max,a,somme);
	  probleme=OUI;
	  break;
	}
    }
  if(probleme==OUI)
    {
      printf("probleme initialisation\n");
    }

  /*calcul de la log-vraisemblance associee a ces parametres*/
  logVraisemblanceCourante=
    calculLogVraisemblanceIndividusSimplex1cluster(nbDimensions,points,cluster,cluster.alphas);

  /*on calcule les parametres par la methode du gradient a pas variable*/
  nbIterations=0;
  while(1)
    {
      normeDeriveeCarre=0;
      
      /*on calcule le gradient*/
      for(j=0;j<nbDimensions;j++)
	{
	  gradient[j]=(digamma(sommeAlphas)-digamma(cluster.alphas[j]))*(double)cluster.nbIndividus;

	  for(k=0;k<cluster.nbIndividus;k++)
	    {
	      gradient[j]+=log(points[k][j]);
	    }
	  normeDeriveeCarre+=gradient[j]*gradient[j];
	}
      
      denominateur=1.0/((double)cluster.nbIndividus*trigamma(sommeAlphas));
      for(j=0;j<nbDimensions;j++)
	{
	  denominateur-=1.0/((double)cluster.nbIndividus*trigamma(cluster.alphas[j]));
	}

      for(j=0;j<nbDimensions;j++)
	{
	  for(k=j;k<nbDimensions;k++)
	    {
	      inverseHessien[j][k]=1.0/(-(double)cluster.nbIndividus*trigamma(cluster.alphas[j])*
					(double)cluster.nbIndividus*trigamma(cluster.alphas[k]))/denominateur;
	      inverseHessien[k][j]=inverseHessien[j][k];
	    }
	  inverseHessien[j][j]-=1.0/((double)cluster.nbIndividus*trigamma(cluster.alphas[j])); 
	}
      
      for(j=0;j<nbDimensions;j++)
	{
	  nouveauxAlphas[j]=cluster.alphas[j];
	  for(k=0;k<nbDimensions;k++)
	    {
	      nouveauxAlphas[j]-=inverseHessien[j][k]*gradient[k];
	    }
	}
      
      /*calcul de la log-vraisemblance associee aux nouveaux parametres*/
      nouvelleLogVraisemblance=
	calculLogVraisemblanceIndividusSimplex1cluster(nbDimensions,points,cluster,nouveauxAlphas);

      differenceAlphas=0;
      for(j=0;j<nbDimensions;j++)
	{
	  differenceAlphas+=(nouveauxAlphas[j]-cluster.alphas[j])*
	    (nouveauxAlphas[j]-cluster.alphas[j]);
	}
      differenceAlphas=sqrt(differenceAlphas);
      
      sommeAlphas=0;
      for(j=0;j<nbDimensions;j++)
	{
	  if(!((cluster.alphas[j]>ZERO_LIMIT)&&(cluster.alphas[j]<999999999.0)))
	    {
	      probleme=OUI;
	      break;
	    }
	} 
      if(probleme==OUI)
	{
	  break;
	} else {
	for(j=0;j<nbDimensions;j++)
	  {
	    cluster.alphas[j]=nouveauxAlphas[j];
	    sommeAlphas+=cluster.alphas[j];
	  }
      }
      
      if(fabs(logVraisemblanceCourante-nouvelleLogVraisemblance)<ZERO_LIMIT)
	{
	  break;
	}
      logVraisemblanceCourante=nouvelleLogVraisemblance;
      nbIterations++;
      if(nbIterations>nbIterationsMax)
	{	  
	  break;
	}
    }
  probleme=NON;
  	
  /*desallocation memoire*/
  ckfree(variances);
  ckfree(gradient);
  ckfree(nouveauxAlphas);
  for(i=0;i<nbDimensions;i++)
    {
      ckfree(inverseHessien[i]);
    }
  ckfree(inverseHessien);
  ckfree(fj);
  /*fin desallocation memoire*/

  if(probleme==OUI)
    {
      printf("probleme\n");
      return OUI;
    }
  else
    {
      return NON;
    }
   
} 

/************************************************************************************************/
/*                                                                                              */
/*Procedure de calcul de la log-vraisemblance d'une classe d'individus pour une loi de Dirichlet*/
/*                                                                                              */
/************************************************************************************************/

double calculLogVraisemblanceIndividusSimplex1cluster(int nbDimensions,double **points,
						      cluster_t cluster,double *alphas)
{
  /*declaration de variables*/
  int i,j;
  double somme,logVraisemblance;
  /*fin declaration de variables*/

  logVraisemblance=0;
  somme=0;
  for(i=0;i<nbDimensions;i++)
    {
      somme+=alphas[i];
    }
  logVraisemblance+=lgammafn(somme);
  for(i=0;i<nbDimensions;i++)
    {
      logVraisemblance-=lgammafn(alphas[i]);
    }
  logVraisemblance*=(double)cluster.nbIndividus;

  for(j=0;j<nbDimensions;j++)
    {
      somme=0;
      for(i=0;i<cluster.nbIndividus;i++)
	{
	  somme+=log(points[i][j]);
	}
      logVraisemblance+=somme*(alphas[j]-1.0);
    }

  return logVraisemblance;
}

/********************************************************************/
/*                                                                  */
/*Procedure de calcul des log-vraisemblances pour tous les individus*/
/*par rapport a tous les clusters                                   */
/*                                                                  */
/********************************************************************/

void calculVraisemblancesPointsSimplex(int nbIndividus,int nbDimensions,individu_t **individus,
					  int nbClusters,cluster_t *clusters) {
  /*declaration de variables*/
  int i,j,k;
  double *logVraisemblances,somme,commun;
  /*fin declaration de variables*/

  /*allocation memoire*/
  logVraisemblances=(double *)ckalloc(sizeof(double)*nbClusters);
  /*fin allocation memoire*/

  for(i=0;i<nbClusters;i++)
    {  
      commun=0;
      somme=0;
      for(j=0;j<nbDimensions;j++)
	{
	  somme+=clusters[i].alphas[j];  
	  commun-=lgammafn(clusters[i].alphas[j]);
	}
      commun+=log((double)clusters[i].nbIndividus/(double)nbIndividus)+lgammafn(somme);
      for(j=0;j<nbIndividus;j++)
	{
	  individus[j]->vrais[i]=commun;
	  for(k=0;k<nbDimensions;k++)
	    {
	      individus[j]->vrais[i]+=(clusters[i].alphas[k]-1.0)*log(individus[j]->valeursTraitees[k]);
	    }
	}
    }

  /*desallocation memoire*/
  ckfree(logVraisemblances);
  /*fin desallocation memoire*/
}



