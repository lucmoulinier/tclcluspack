#include "mainClusPack.h"
#include "divide.h"
#include "tools.h"
#include "2means.h"
#include <tcl.h>

/******************************************************/
/*                                                    */
/*Procedure de realisation d'un clustering d'un groupe*/ 
/*en deux par la methode des k-means avec k=2         */
/*                                                    */
/******************************************************/

int twoMeans(int nbDimensions,cluster_t *groupe,cluster_t *groupe1,cluster_t *groupe2,
	      int typeDonnees)
{
  /*declaration des variables*/
  int i,j,point1,point2,*groupesDesIndividus,succesAffectation;
  int compteurIterations,nbIterationsTwoMeans,compteurWhile;
  double *centre,*centre1,*centre2,inertieInterClassePrecedente;
  double inertieInterClasseCourante,meilleureInertieInterClasse;
  /*fin declaration des variables*/

  nbIterationsTwoMeans=10;
  meilleureInertieInterClasse=NOT_DEFINED;
  
  /*allocation memoire*/
  centre1=(double *)ckalloc(sizeof(double)*nbDimensions);
  centre2=(double *)ckalloc(sizeof(double)*nbDimensions);
  centre=(double *)ckalloc(sizeof(double)*nbDimensions);
  groupesDesIndividus=(int *)ckalloc(sizeof(int)*groupe->nbIndividus);
  /*fin allocation memoire*/

  /*  if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
    {
      for(i=0;i<groupe->nbIndividus;i++)
	{
	  centre[i]=0;
	  for(j=0;j<groupe->nbIndividus;j++)
	    {
	      centre[i]+=groupe->individus[j]->valeursTraitees[groupe->individus[i]->id];
	    }
	  centre[i]/=(double)groupe->nbIndividus;
	}
    }
  else
  {*/
      for(i=0;i<nbDimensions;i++)
	{
	  centre[i]=0;
	  for(j=0;j<groupe->nbIndividus;j++)
	    {
	      centre[i]+=groupe->individus[j]->valeursTraitees[i];
	    }
	  centre[i]/=(double)groupe->nbIndividus;
	}
      /*}*/

  for(compteurIterations=0;compteurIterations<nbIterationsTwoMeans;compteurIterations++)
    {
      point1=(int)floor((double)groupe->nbIndividus-1.0)*runif(0,1);
      compteurWhile=0;
      while(compteurWhile<50)
	{
	  point2=(int)floor((double)groupe->nbIndividus-1.0)*runif(0,1);
	  if(calculDistance(nbDimensions,groupe->individus[point1]->valeursTraitees,
			     groupe->individus[point2]->valeursTraitees)>ZERO_LIMIT)
	    {
	      break;
	    }
	  compteurWhile++;
	}
      if(compteurWhile==50)
	{
	  ckfree((char *)centre1);
	  ckfree((char *)centre2);
	  ckfree((char *)centre);
	  ckfree((char *)groupesDesIndividus);

	  return NON;
	}

      /*initialisation*/
      /*      if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
	{
	  for(i=0;i<groupe->nbIndividus;i++)
	    {
	      centre1[i]=
		groupe->individus[point1]->valeursTraitees[groupe->individus[i]->id];
	      centre2[i]=
		groupe->individus[point2]->valeursTraitees[groupe->individus[i]->id];
	    }
	}
      else
      {*/
	  for(i=0;i<nbDimensions;i++)
	    {
	      centre1[i]=groupe->individus[point1]->valeursTraitees[i];
	      centre2[i]=groupe->individus[point2]->valeursTraitees[i];
	    }
	  /*}*/

      /*premiere affectation des individus aux deux groupes*/
      succesAffectation=
	affectationIndividus2means(nbDimensions,groupe->nbIndividus,groupe->individus,
				    centre1,centre2,groupesDesIndividus,
				    &inertieInterClassePrecedente,typeDonnees,centre);
      if(succesAffectation==NON)
	{
	  return NON;
	}

      /*fin initialisation*/
      compteurWhile=0;
      while(compteurWhile<10)
	{
	  succesAffectation=
	    affectationIndividus2means(nbDimensions,groupe->nbIndividus,groupe->individus,
					centre1,centre2,groupesDesIndividus,
					&inertieInterClasseCourante,typeDonnees,centre);
	  if(succesAffectation==NON)
	    {
	      return NON;
	    }

	  if((inertieInterClasseCourante-inertieInterClassePrecedente<ZERO_LIMIT)&&
	     (-inertieInterClasseCourante+inertieInterClassePrecedente<ZERO_LIMIT))
	    {
	      break;
	    }
	  else
	    {
	      inertieInterClassePrecedente=inertieInterClasseCourante;
	    }
	  compteurWhile++;
	}
      
      if((meilleureInertieInterClasse==NOT_DEFINED)||
	 (inertieInterClasseCourante>meilleureInertieInterClasse))
	{
	  groupe1->nbIndividus=0;
	  groupe2->nbIndividus=0;
	  for(i=0;i<groupe->nbIndividus;i++)
	    {
	      if(groupesDesIndividus[i]==GROUPE1)
		{
		  groupe1->individus[groupe1->nbIndividus]=groupe->individus[i];
		  (groupe1->nbIndividus)++;
		}
	      else
		{
		  groupe2->individus[groupe2->nbIndividus]=groupe->individus[i];
		  (groupe2->nbIndividus)++;
		}
	    }
	  meilleureInertieInterClasse=inertieInterClasseCourante;
	}
    }      

  /*desallocation memoire*/
  ckfree((char *)centre1);
  ckfree((char *)centre2);
  ckfree((char *)centre);
  ckfree((char *)groupesDesIndividus);
  /*fin desallocation memoire*/

  return OUI;
}


/***************************************************/
/*                                                 */
/*Procedure d'affectation des individus aux groupes*/
/*et de calcul des nouveaux centres de gravite     */
/*                                                 */
/***************************************************/

int affectationIndividus2means(int nbDimensions,int nbIndividus,individu_t **individus,
				double *centre1,double *centre2,int *groupesDesIndividus,
				double *inertieInterClasse,int typeDonnees,double *centre)
{
  /*declaration des variables*/
  int i,j,nbIndividusGroupe1,nbIndividusGroupe2;
  double distance1,distance2,contribution1,contribution2;
  /*fin declaration des variables*/

  nbIndividusGroupe1=0;
  nbIndividusGroupe2=0;
  for(i=0;i<nbIndividus;i++)
    {
      /*      if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
	{
	  distance1=0;
	  distance2=0;
	  for(j=0;j<nbIndividus;j++)
	    {
	      distance1+=(centre1[j]-individus[i]->valeursTraitees[individus[j]->id])*
		(centre1[j]-individus[i]->valeursTraitees[individus[j]->id]);
	      distance2+=(centre2[j]-individus[i]->valeursTraitees[individus[j]->id])*
		(centre2[j]-individus[i]->valeursTraitees[individus[j]->id]);
	    }
	}
      else
      {*/
	  distance1=calculDistanceCarre(nbDimensions,centre1,individus[i]->valeursTraitees);
	  distance2=calculDistanceCarre(nbDimensions,centre2,individus[i]->valeursTraitees);
	  /*}*/

      if(distance1<distance2)
	{
	  groupesDesIndividus[i]=GROUPE1;
	  nbIndividusGroupe1++;
	}
      else
	{
	  groupesDesIndividus[i]=GROUPE2;
	  nbIndividusGroupe2++;
	}
    }
  if((nbIndividusGroupe1==0)||(nbIndividusGroupe2==0))
    {
      return NON;
    }

  for(i=0;i<nbDimensions;i++)
    {
      centre1[i]=0;
      centre2[i]=0;
    }
  
  for(i=0;i<nbIndividus;i++)
    {
      if(groupesDesIndividus[i]==GROUPE1)
	{
	  /*	  if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
	    {
	      for(j=0;j<nbIndividus;j++)
		{
		  centre1[j]+=individus[i]->valeursTraitees[individus[j]->id];
		}
	    }
	  else
	  {*/
	      for(j=0;j<nbDimensions;j++)
		{
		  centre1[j]+=individus[i]->valeursTraitees[j];
		}
	      /*}*/
	}
      else
	{
	  /*	  if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
	    {
	      for(j=0;j<nbIndividus;j++)
		{
		  centre2[j]+=individus[i]->valeursTraitees[individus[j]->id];
		}
	    }
	  else
	  {*/
	      for(j=0;j<nbDimensions;j++)
		{
		  centre2[j]+=individus[i]->valeursTraitees[j];
		}
	      /*}*/
	}
    }
 
  /*  if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
    {
      for(i=0;i<nbIndividus;i++)
	{
	  centre1[i]/=(double)nbIndividusGroupe1;
	  centre2[i]/=(double)nbIndividusGroupe2;
	}
    }
  else
  {*/
      for(i=0;i<nbDimensions;i++)
	{
	  centre1[i]/=(double)nbIndividusGroupe1;
	  centre2[i]/=(double)nbIndividusGroupe2;
	}
      /* }

  if((typeDonnees==SIMILARITES)||(typeDonnees==ALIGNEMENT))
    { 
      contribution1=0;
      for(j=0;j<nbIndividus;j++)
	{
	  contribution1+=(centre[j]-centre1[j])*(centre[j]-centre1[j]);
	}
      contribution2=0;
      for(j=0;j<nbIndividus;j++)
	{
	  contribution2+=(centre[j]-centre2[j])*(centre[j]-centre2[j]);
	}
    }
  else
  {*/
      contribution1=0;
      for(j=0;j<nbDimensions;j++)
	{
	  contribution1+=(centre[j]-centre1[j])*(centre[j]-centre1[j]);
	}
      contribution2=0;
      for(j=0;j<nbDimensions;j++)
	{
	  contribution2+=(centre[j]-centre2[j])*(centre[j]-centre2[j]);
	}
      /*}*/
  *inertieInterClasse=contribution1*nbIndividusGroupe1+contribution2*nbIndividusGroupe2;

  return OUI;
}

