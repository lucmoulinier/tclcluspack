#include "mainClusPack.h"
#include "tools.h"
#include "dissimilarities_clustering.h"
#include "tcl.h"

/*************************************************************/
/*                                                           */
/*classification des dissimilarites pour trouver le nombre de*/
/*dissimilarites elevees et donc le nombre de groupes        */
/*(avec l'exposant de Holder)                                */
/*                                                           */
/*************************************************************/

void clusteringDissimilaritesHolder(int nbDissimilarites,
				    double *historiqueDissimilarites,
				    int *nbClusters,
				    int resolution)
{
  /*declaration de variables*/
  int i;
  int auMoinsUnExposantHolderUnitaire,nbExposantsHolder,nbCoupuresCandidates;
  int *coupuresCandidates,coupureSelectionnees,coupureParDefaut;
  int nbDissimilaritesElevees=0;
  double *exposantsHolder,*exposantsHolderOrdonnes;
  /*fin declaration de variables*/

  /*allocation memoire*/
  exposantsHolder=(double *)ckalloc(sizeof(double)*nbDissimilarites);
  exposantsHolderOrdonnes=(double *)ckalloc(sizeof(double)*nbDissimilarites);
  coupuresCandidates=(int *)ckalloc(sizeof(int)*nbDissimilarites);
  /*fin allocation memoire*/

  calculExposantsHolder(nbDissimilarites,
			historiqueDissimilarites,
			&nbExposantsHolder,
			exposantsHolder);
  
  auMoinsUnExposantHolderUnitaire=NON;
  for(i=0;i<nbExposantsHolder;i++)
    {
      if(exposantsHolder[i]>EXPOSANT_HOLDER_SEUIL)
	{
	  auMoinsUnExposantHolderUnitaire=OUI;
	  break;
	}
    }
  
  if(auMoinsUnExposantHolderUnitaire==OUI)
    {
      /*calcul du nombre de divisions pour la dimension fractale*/
      
      /*clustering proprement dit des dissimilarites*/
      separationDissimilaritesEvolue(nbExposantsHolder,
				     nbDissimilarites,
				     historiqueDissimilarites,
				     &nbDissimilaritesElevees,
				     exposantsHolder);      
      
      nbCoupuresCandidates=0;
      coupuresCandidates[0]=0;
      
      if(resolution!=0)
	{
	  for(i=0;i<nbExposantsHolder-1;i++)
	    {
	      if((exposantsHolder[i]<=EXPOSANT_HOLDER_SEUIL)&&
		 (exposantsHolder[i+1]>EXPOSANT_HOLDER_SEUIL))
		{
		  coupuresCandidates[nbCoupuresCandidates]=i;
		  if((i+1)==nbDissimilaritesElevees)
		    {
		      coupureParDefaut=nbCoupuresCandidates;
		    }
		  nbCoupuresCandidates++;
		}
	    }

	  printf("coupure_par_defaut : %d\n",coupureParDefaut);

	  coupureSelectionnees=resolution+coupureParDefaut;
	  if(coupureSelectionnees<0)
	    {
	      nbDissimilaritesElevees=coupuresCandidates[0]+1;
	    }
	  else if(coupureSelectionnees>=nbCoupuresCandidates)
	    {
	      nbDissimilaritesElevees=coupuresCandidates[nbCoupuresCandidates-1]+1;
	    }
	  else
	    {
	      nbDissimilaritesElevees=coupuresCandidates[coupureSelectionnees]+1;
	    }
	}
    }
  else
    {
      /*on separe les valeurs elevees des valeurs faibles*/
      separationDissimilaritesBete(nbDissimilarites,
				   historiqueDissimilarites,
				   &nbDissimilaritesElevees);
    }
  
  *nbClusters=nbDissimilaritesElevees+1;

  /*desallocation memoire*/
  ckfree((char *)exposantsHolder);
  ckfree((char *)exposantsHolderOrdonnes);
  ckfree((char *)coupuresCandidates);
  /*fin desallocation memoire*/

  return;
}

/******************************************************************************/
/*                                                                            */
/*Procedure de calcul des exposants d'Holder sur les valeurs de dissimilarites*/
/*                                                                            */
/******************************************************************************/

void calculExposantsHolder(int nbDissimilarites,double *historiqueDissimilarites,
			   int *nbExposantsHolder,double *exposantsHolder)
{
  /*declaration des variables*/
  int i,j,nbMesures,compteurExposantsHolder;
  double **valeursPourRegression,mesureMax;
  /*fin declaration des variables*/

  /*allocation memoire*/
  valeursPourRegression=(double **)ckalloc(sizeof(double *)*(nbDissimilarites-1));
  for(i=0;i<nbDissimilarites-1;i++)
    {
      valeursPourRegression[i]=(double *)ckalloc(sizeof(double)*2);
    }
  /*fin allocation memoire*/

  mesureMax=historiqueDissimilarites[0]-historiqueDissimilarites[nbDissimilarites-1];
  
  compteurExposantsHolder=0;
  for(i=0;i<nbDissimilarites-1;i++)
    {
      nbMesures=(int)floor(log((double)(nbDissimilarites-i))/log(2.0));
      for(j=nbMesures;j>0;j--)
	{
	  if(((historiqueDissimilarites[i]-
	      historiqueDissimilarites[i+(int)pow(2.0,(double)j)-1])/mesureMax)<ZERO_LIMIT)
	    {
	      nbMesures--;
	    }
	  else
	    {
	      valeursPourRegression[nbMesures-j][1]=
		log((double)(historiqueDissimilarites[i]-
		     historiqueDissimilarites[i+(int)pow(2.0,(double)j)-1])/mesureMax);
	      valeursPourRegression[nbMesures-j][0]=-log(2.0)*(double)(nbMesures-j);
	    }
	}

      if(nbMesures>1)
	{
	  exposantsHolder[compteurExposantsHolder]=
	    regressionLineaire(nbMesures,valeursPourRegression);
	  if(exposantsHolder[compteurExposantsHolder]>1.0)
	    {
	      exposantsHolder[compteurExposantsHolder]=1.0;
	    }
	  compteurExposantsHolder++;
	}
    }

  *nbExposantsHolder=compteurExposantsHolder;

  /*desallocation memoire*/
  for(i=0;i<nbDissimilarites-1;i++)
    {
      ckfree((char *)valeursPourRegression[i]);
    }
  ckfree(valeursPourRegression);
  /*fin desallocation memoire*/
}

/**********************************************************************/
/*                                                                    */
/*Procedure de separation en deux groupes des valeurs de dissimilarite*/
/*par calcul du minimum d'inertie intraclasse                         */
/*                                                                    */
/**********************************************************************/

void separationDissimilaritesBete(int nbDissimilarites,double *historiqueDissimilarites,
				    int *nbDissimilaritesPremierGroupe)
{
  /*declaration des variables*/
  int i,j,nbDissimilaritesValeurMinimum,nbDissimilarites1,nbDissimilarites2;
  double valeurCourante,valeurMinimum,centreGravite1,centreGravite2;
  /*fin declaration des variables*/

  valeurMinimum=-1;
  nbDissimilaritesValeurMinimum=1;
  for(i=1;i<nbDissimilarites;i++)
    {
      nbDissimilarites1=i;
      nbDissimilarites2=nbDissimilarites-nbDissimilarites1;
      
      centreGravite1=0;
      for(j=0;j<nbDissimilarites1;j++)
	{
	  centreGravite1+=historiqueDissimilarites[j];
	}
      centreGravite1/=(double)nbDissimilarites1;
    
      centreGravite2=0;
      for(j=nbDissimilarites1;j<nbDissimilarites;j++)
	{
	  centreGravite2+=historiqueDissimilarites[j];
	}
      centreGravite2/=(double)nbDissimilarites2;

      valeurCourante=0;

      for(j=0;j<nbDissimilarites1;j++)
	{
	  valeurCourante+=((double)(historiqueDissimilarites[j]-centreGravite1))*
	    ((double)(historiqueDissimilarites[j]-centreGravite1));
	}
      for(j=nbDissimilarites1;j<nbDissimilarites;j++)
	{
	  valeurCourante+=((double)(historiqueDissimilarites[j]-centreGravite2))*
	    ((double)(historiqueDissimilarites[j]-centreGravite2));
	}

      if((valeurMinimum<0)||(valeurCourante<valeurMinimum))
	{
	  valeurMinimum=valeurCourante;
	  nbDissimilaritesValeurMinimum=nbDissimilarites1;
	}
    }

  *nbDissimilaritesPremierGroupe=nbDissimilaritesValeurMinimum;
}

/************************************************************************/
/*                                                                      */
/*Procedure de separation en deux groupes des valeurs de dissimilarite  */
/*par calcul du minimum d'inertie intraclasse mais uniquement aux points*/
/*d'exposants d'Holder faibles                                          */
/*                                                                      */
/************************************************************************/

void separationDissimilaritesEvolue(int nbExposantsHolder,int nbDissimilarites,
				      double *historiqueDissimilarites,
				      int *nbDissimilaritesPremierGroupe,
				      double *exposantsHolder)
{
  /*declaration des variables*/
  int i,j,nbDissimilaritesValeurMinimum,nbDissimilarites1,nbDissimilarites2,compteur;
  double valeurCourante,valeurMinimum,centreGravite1,centreGravite2;
  double *inertiesIntraClasse;
  /*fin declaration des variables*/

  /*allocation memoire*/
  inertiesIntraClasse=(double *)ckalloc(sizeof(double)*nbDissimilarites);
  /*fin allocation memoire*/

  valeurMinimum=-1;
  nbDissimilaritesValeurMinimum=1;
  compteur=0;
  for(i=1;i<nbExposantsHolder;i++)
    {
      if(i>nbDissimilarites/2)
	{
	  break;
	}
      if((exposantsHolder[i-1]<=EXPOSANT_HOLDER_SEUIL)&&
	 (exposantsHolder[i]>EXPOSANT_HOLDER_SEUIL))
	{
	  nbDissimilarites1=i;
	  nbDissimilarites2=nbDissimilarites-nbDissimilarites1;
	  
	  centreGravite1=0;
	  for(j=0;j<nbDissimilarites1;j++)
	    {
	      centreGravite1+=historiqueDissimilarites[j];
	    }
	  centreGravite1/=(double)nbDissimilarites1;
	  
	  centreGravite2=0;
	  for(j=nbDissimilarites1;j<nbDissimilarites;j++)
	    {
	      centreGravite2+=historiqueDissimilarites[j];
	    }
	  centreGravite2/=(double)nbDissimilarites2;
	  
	  valeurCourante=0;
	  
	  for(j=0;j<nbDissimilarites1;j++)
	    {
	      valeurCourante+=((double)(historiqueDissimilarites[j]-centreGravite1))*
		((double)(historiqueDissimilarites[j]-centreGravite1));
	    }
	  for(j=nbDissimilarites1;j<nbDissimilarites;j++)
	    {
	      valeurCourante+=((double)(historiqueDissimilarites[j]-centreGravite2))*
		((double)(historiqueDissimilarites[j]-centreGravite2));
	    }
	  
	  if((valeurMinimum<0)||(valeurCourante<valeurMinimum))
	    {
	      valeurMinimum=valeurCourante;
	      nbDissimilaritesValeurMinimum=nbDissimilarites1;
	    }
	  /*printf("inertie intra classe en %d : %.2f\n",i,valeurCourante);*/
	  inertiesIntraClasse[i]=valeurCourante;
	  compteur++;
	}
    }

  *nbDissimilaritesPremierGroupe=nbDissimilaritesValeurMinimum;

  /*desallocation memoire*/
  ckfree((char *)inertiesIntraClasse);
  /*fin desallocation memoire*/
}

// 06 23 84 31 88
