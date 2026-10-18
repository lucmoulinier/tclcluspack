#include "mainClusPack.h"
#include "divide.h"
#include "space.h"
#include "tools.h"

/**************************************************************/
/*                                                            */
/*Procedure de test pour savoir si l'on doit scinder un groupe*/
/*par la methode DPC (Density of Points Clustering)           */ 
/*                                                            */
/**************************************************************/

int testDPC(int nbIndividus,int nbDimensions,cluster_t *groupePere,
	    cluster_t *groupe1,cluster_t *groupe2,double **coordonneesIndividus,
	    double *cotesHyperpave,double *valeursMinimum,double risque,
	    int typeDonnees,int typeDensite,int donneesNormalisees,int nbSimulations)
{
  /*declaration de variables*/
  int i,j,k,individu1,individu2,nbIterations,nbIterations2;
  int nbFoisHypotheseRejetee,compteurSimulations,compteurPopulationReference,compteur;
  double *individuFictif,exposantHolder,*exposantsHolderPere;
  double moyennePere,moyenneAleatoire,moyenne;
  double *exposantsHolderAleatoire,varianceResiduelle,varianceFacteur,distanceCourante;
  double distanceMinimum,f,P,*directionAleatoire,signe,m,s,varianceResiduellePere;
  /*fin declaration de variables*/

  nbFoisHypotheseRejetee=0;

  nbIterations2=groupe1->nbIndividus*groupe2->nbIndividus/2;

  if(nbIterations2>500)
    {
      nbIterations2=500;
    }
  else if(nbIterations2<2)
    {
      nbIterations2=2;
    }

  /*allocation memoire*/
  individuFictif=(double *)ckalloc(sizeof(double)*nbDimensions);
  directionAleatoire=(double *)ckalloc(sizeof(double)*nbDimensions);
  exposantsHolderPere=(double *)ckalloc(sizeof(double)*(2*groupePere->nbIndividus+2));
  exposantsHolderAleatoire=(double *)ckalloc(sizeof(double)*nbIterations2);
  /*fin allocation memoire*/

  /*printf("taille1 : %d ; taille2 : %d\n",groupe1->nbIndividus,groupe2->nbIndividus);*/
  /*for(i=0;i<nbDimensions;i++)
    {
    individuFictif[i]=0.45;
    }
    
    calculExposantsHolderPoints(nbDimensions,nbIndividus,coordonneesIndividus,
    1,&individuFictif,&exposantHolder,
    cotesHyperpave,valeursMinimum,NON);
    for(i=0;i<nbIndividus;i++)
    {
    for(j=0;j<nbDimensions;j++)
    {
    printf("%.3f ",coordonneesIndividus[i][j]);
    }
    printf("\n");
    }
    printf("exposant holder : %f\n",exposantHolder);
    exit(0);*/

  moyennePere=0;
  compteurPopulationReference=0;
  
  if(groupe1->nbIndividus>1)
    {
      if(groupe1->nbIndividus<50)
	{
	  nbIterations=groupe1->nbIndividus;
	}
      else
	{
	  nbIterations=50;
	}
      for(i=0;i<nbIterations;i++)
	{
	  individu1=(int)floor((double)groupe1->nbIndividus*runif(0,1));
	  if(individu1==groupe1->nbIndividus)
	    {
	      individu1=groupe1->nbIndividus-1;
	    }

	  distanceMinimum=-1.0;
	  for(j=0;j<groupe1->nbIndividus;j++)
	    {
	      if(j!=individu1)
		{
		  distanceCourante=0;
		  for(k=0;k<nbDimensions;k++)
		    {
		      distanceCourante+=
			(groupe1->individus[individu1]->valeursTraitees[k]-
			 groupe1->individus[j]->valeursTraitees[k])*
			(groupe1->individus[individu1]->valeursTraitees[k]-
			 groupe1->individus[j]->valeursTraitees[k]);
		    }
		  if((distanceMinimum<0)||(distanceCourante<distanceMinimum))
		    {
		      distanceMinimum=distanceCourante;
		      individu2=j;
		    }
		}
	    }

	  for(j=0;j<nbDimensions;j++)
	    {
	      directionAleatoire[j]=groupe1->individus[individu2]->valeursTraitees[j]-
		groupe1->individus[individu1]->valeursTraitees[j];
	    }
	 
	  compteur=0;
	  for(j=0;j<nbDimensions;j++) 
	    { 
	      if(runif(0,1)<0.5)
		{
		  signe=-1.0;
		  compteur++;
		}
	      else
		{
		  signe=1.0;
		}
	      /*if((compteur>1)&&(compteur>nbDimensions/5)&&(typeDensite==DENSITE2))
		{
		  signe=1.0;
		  }*/
	      individuFictif[j]=groupe1->individus[individu1]->valeursTraitees[j]+
		signe*directionAleatoire[j];
	    }

	  if(donneesNormalisees==OUI)
	    {
	      /*on normalise l'individu fictif*/
	      m=0;
	      for(j=0;j<nbDimensions;j++)
		{
		  m+=individuFictif[j];
		}
	      m/=(double)nbDimensions;

	      s=0;
	      for(j=0;j<nbDimensions;j++)
		{
		  s+=((double)(m-individuFictif[j]))*((double)(m-individuFictif[j]));
		}
	      s/=(double)nbDimensions;
	      s=sqrt((double)s);

	      for(j=0;j<nbDimensions;j++)
		{
		  individuFictif[j]=(individuFictif[j]-m)/s;
		}
	    }

	  /*on verifie que l'individu fictif n'est pas hors-norme*/
	  for(j=0;j<nbDimensions;j++)
	    {
	      if(individuFictif[j]<valeursMinimum[j])
		{
		  individuFictif[j]=valeursMinimum[j];
		}
	      else if(individuFictif[j]>valeursMinimum[j]+cotesHyperpave[j])
		{
		  individuFictif[j]=valeursMinimum[j]+cotesHyperpave[j];
		}
	    }
	  
	  calculExposantsHolderPoints(nbDimensions,nbIndividus,coordonneesIndividus,
					 1,&individuFictif,&exposantHolder,
					 cotesHyperpave,valeursMinimum,NON);

	  exposantsHolderPere[compteurPopulationReference]=exposantHolder;
	  moyennePere+=exposantHolder;
	  compteurPopulationReference++;
	}
    }

  if(groupe2->nbIndividus>1)
    {
      if(groupe2->nbIndividus<50)
	{
	  nbIterations=groupe2->nbIndividus;
	}
      else
	{
	  nbIterations=50;
	}
      for(i=0;i<nbIterations;i++)
	{
	  individu1=(int)floor((double)groupe2->nbIndividus*runif(0,1));
	  if(individu1==groupe2->nbIndividus)
	    {
	      individu1=groupe2->nbIndividus-1;
	    }  
	  distanceMinimum=-1.0;
	  for(j=0;j<groupe2->nbIndividus;j++)
	    {
	      if(j!=individu1)
		{
		  distanceCourante=0;
		  for(k=0;k<nbDimensions;k++)
		    {
		      distanceCourante+=
			(groupe2->individus[individu1]->valeursTraitees[k]-
			 groupe2->individus[j]->valeursTraitees[k])*
			(groupe2->individus[individu1]->valeursTraitees[k]-
			 groupe2->individus[j]->valeursTraitees[k]);
		    }
		  if((distanceMinimum<0)||(distanceCourante<distanceMinimum))
		    {
		      distanceMinimum=distanceCourante;
		      individu2=j;
		    }
		}
	    }

	  for(j=0;j<nbDimensions;j++)
	    {
	      directionAleatoire[j]=groupe2->individus[individu2]->valeursTraitees[j]-
		groupe2->individus[individu1]->valeursTraitees[j];
	    }

	  compteur=0;
	  for(j=0;j<nbDimensions;j++) 
	    { 
	      if(runif(0,1)<0.5)
		{
		  signe=-1.0;
		  compteur++;
		}
	      else
		{
		  signe=1.0;
		}
	      /*if((compteur>1)&&(compteur>nbDimensions/5)&&(typeDensite==DENSITE2))
		{
		  signe=1.0;
		  }*/
	      individuFictif[j]=groupe2->individus[individu1]->valeursTraitees[j]+
		signe*directionAleatoire[j];
	    }

	  if(donneesNormalisees==OUI)
	    {
	      /*on normalise l'individu fictif*/
	      m=0;
	      for(j=0;j<nbDimensions;j++)
		{
		  m+=individuFictif[j];
		}
	      m/=(double)nbDimensions;

	      s=0;
	      for(j=0;j<nbDimensions;j++)
		{
		  s+=((double)(m-individuFictif[j]))*((double)(m-individuFictif[j]));
		}
	      s/=(double)nbDimensions;
	      s=sqrt((double)s);

	      for(j=0;j<nbDimensions;j++)
		{
		  individuFictif[j]=(individuFictif[j]-m)/s;
		}
	    }

	  /*on verifie que l'individu fictif n'est pas hors-norme*/
	  for(j=0;j<nbDimensions;j++)
	    {
	      if(individuFictif[j]<valeursMinimum[j])
		{
		  individuFictif[j]=valeursMinimum[j];
		}
	      else if(individuFictif[j]>valeursMinimum[j]+cotesHyperpave[j])
		{
		  individuFictif[j]=valeursMinimum[j]+cotesHyperpave[j];
		}
	    }
	  	  
	  calculExposantsHolderPoints(nbDimensions,nbIndividus,coordonneesIndividus,
					 1,&individuFictif,&exposantHolder,
					 cotesHyperpave,valeursMinimum,NON);

	  exposantsHolderPere[compteurPopulationReference]=exposantHolder;
	  moyennePere+=exposantHolder;
	  compteurPopulationReference++;
	}
    }

  moyennePere/=(double)(compteurPopulationReference);
  varianceResiduellePere=0;
  for(i=0;i<compteurPopulationReference;i++)
    { 
      varianceResiduellePere+=((double)(exposantsHolderPere[i]-moyennePere))*
	((double)(exposantsHolderPere[i]-moyennePere));
    } 

  for(compteurSimulations=0;compteurSimulations<nbSimulations;compteurSimulations++)
    {
      moyenneAleatoire=0;
      for(i=0;i<nbIterations2;i++)
	{
	  /*selection aleatoire de deux individus du groupe pere*/
	  individu1=(int)floor((double)(groupe1->nbIndividus-1)*runif(0,1));
	  individu2=(int)floor((double)(groupe2->nbIndividus-1)*runif(0,1));
	
	  for(j=0;j<nbDimensions;j++) 
	    { 
	      if((typeDonnees==ALIGNEMENT)||(typeDensite==DENSITE1))
		{
		  individuFictif[j]=(groupe1->individus[individu1]->valeursTraitees[j]-
				     groupe2->individus[individu2]->valeursTraitees[j])*0.5
		    +groupe2->individus[individu2]->valeursTraitees[j];
		}
	      else
		{
		  individuFictif[j]=(groupe1->individus[individu1]->valeursTraitees[j]-
				     groupe2->individus[individu2]->valeursTraitees[j])*runif(0,1)
		    +groupe2->individus[individu2]->valeursTraitees[j];
		}
	    }
	  
	  if(donneesNormalisees==OUI)
	    {
	      /*on normalise l'individu fictif*/
	      m=0;
	      for(j=0;j<nbDimensions;j++)
		{
		  m+=individuFictif[j];
		}
	      m/=(double)nbDimensions;
	      
	      s=0;
	      for(j=0;j<nbDimensions;j++)
		{
		  s+=((double)(m-individuFictif[j]))*((double)(m-individuFictif[j]));
		}
	      s/=(double)nbDimensions;
	      s=sqrt((double)s);
	      
	      for(j=0;j<nbDimensions;j++)
		{
		  individuFictif[j]=(individuFictif[j]-m)/s;
		}
	    }

	  calculExposantsHolderPoints(nbDimensions,nbIndividus,coordonneesIndividus,
					 1,&individuFictif,&exposantHolder,
					 cotesHyperpave,valeursMinimum,NON);
	  
	  exposantsHolderAleatoire[i]=exposantHolder;
	  moyenneAleatoire+=exposantHolder;
	}
     
      moyenneAleatoire=moyenneAleatoire/(double)nbIterations2;

      /*realisation du test*/
      moyenne=(moyenneAleatoire*(double)nbIterations2+moyennePere*
	       (double)compteurPopulationReference)
	/(double)(nbIterations2+compteurPopulationReference);
    
      varianceFacteur=(((double)(moyennePere-moyenne))*((double)(moyennePere-moyenne))*
			(double)compteurPopulationReference+((double)(moyenneAleatoire-moyenne))*
			((double)(moyenneAleatoire-moyenne))*(double)nbIterations2)/
	(double)(nbIterations2+compteurPopulationReference);
     
      varianceResiduelle=varianceResiduellePere;
      for(i=0;i<nbIterations2;i++)
	{
	  varianceResiduelle+=((double)(exposantsHolderAleatoire[i]-moyenneAleatoire))*
	    ((double)(exposantsHolderAleatoire[i]-moyenneAleatoire));
	}

      varianceResiduelle/=(double)(nbIterations2+compteurPopulationReference);
      

      if(varianceResiduelle>ZERO_LIMIT)
	{
	  f=varianceFacteur*(double)(nbIterations2+compteurPopulationReference-2)/
	    varianceResiduelle;
	  
	  P=F_FisherSnedecor(1.0,(double)(nbIterations2+compteurPopulationReference-2),
			      (double)f);

	  /*printf("1-P : %f ; %f ; %f\n",1.0-P,moyennePere,moyenneAleatoire);*/
	    
	  if((1.0-P<risque)&&(moyennePere-ZERO_LIMIT>moyenneAleatoire))
	    {  
	      nbFoisHypotheseRejetee++;
	    } 
	}
    }

  /*desallocation memoire*/
  ckfree(individuFictif);
  ckfree(directionAleatoire);
  ckfree(exposantsHolderPere);
  ckfree(exposantsHolderAleatoire);
  /*fin desallocation memoire*/

  if(nbFoisHypotheseRejetee>nbSimulations/2)
    {
      return OUI;
    }
  else
    {
      return NON;
    }
}

