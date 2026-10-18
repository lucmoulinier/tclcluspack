#include "mainClusPack.h"
#include "matrice.h"
#include "emSphericalGaussian.h"
#include "space.h"

/********************************************/
/*                                          */
/*Procedure de calcul de la densite uniforme*/
/*                                          */
/********************************************/

double calculLogDensiteUniformeGaussienneSpherique(int nbIndividus,int nbDimensions,individu_t **individus) {
  /*declaration de variables*/
  double logDensiteLoiUniforme,surfaceSphere;
  /*fin declaration de variables*/
  
  calculSurfaceSphere(nbDimensions-1,&surfaceSphere);

  logDensiteLoiUniforme=-log(surfaceSphere);

  return logDensiteLoiUniforme;
}

/***********************************************/
/*                                             */
/*Procedure de maximisation de la vraisemblance*/
/*                                             */
/***********************************************/

void maximisationVraisemblanceGaussienneSpherique(int nbPoints,double **points,int nbDimensions,
						  cluster_t cluster)
{
  /*declaration de variables*/
  int i;
  double **pointsRot,**pointsProjetes,*angles;
  /*fin declaration de variables*/

  /*allocation memoire*/
  pointsProjetes=(double **)ckalloc(sizeof(double *)*nbPoints);
  for(i=0;i<nbPoints;i++)
    {
      pointsProjetes[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  pointsRot=(double **)ckalloc(sizeof(double *)*nbPoints);
  for(i=0;i<nbPoints;i++)
    {
      pointsRot[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  angles=(double *)ckalloc(sizeof(double)*(nbDimensions-1));
  /*fin allocation memoire*/

  /*calcul du mu moyen*/
  /*generation d'erreurs ici*/ 
  calculMuMoyenGaussienneSpherique(nbPoints,points,nbDimensions,cluster.mu);
  
  /*recuperation des angles de mu*/
  cartesiennes2polaires(nbDimensions,cluster.mu,angles);

  /*calcul des nouvelles coordonnees apres rotation du repere*/
  calculCoordonneesApresRotation(nbDimensions,nbPoints,points,angles,
				 pointsRot); 
   
  /*calcul des coordonnees des points projetes*/
  calculCoordonneesApresProjection(nbDimensions,nbPoints,pointsRot,
				   pointsProjetes);

  /*estimation de la matrice de covariance*/
  calculMatriceCovariance(nbDimensions-1,nbPoints,pointsProjetes,cluster.sigma);

  /*desallocation memoire*/
  for(i=0;i<nbPoints;i++)
    {
      ckfree(pointsProjetes[i]);
    }
  ckfree(pointsProjetes);
  for(i=0;i<nbPoints;i++)
    { 
      ckfree(pointsRot[i]);
    }
  ckfree(pointsRot);
  ckfree(angles);
  /*fin desallocation memoire*/
} 

/**************************************************/
/*                                                */
/*Procedure de calcul du mu moyen pour initialiser*/
/*la maximisation de la vraisemblance             */
/*                                                */
/**************************************************/

void calculMuMoyenGaussienneSpherique(int nbPoints,double **coordonnees,int nbDimensions,
				      double *mu)
{
  /*declaration de variables*/
  int i,j;
  double norme;
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

  norme=0;
  for(j=0;j<nbDimensions;j++)
    {
      mu[j]/=(double)nbPoints;
      norme+=mu[j]*mu[j];
    }
  norme=sqrt(norme);

  for(j=0;j<nbDimensions;j++)
    {
      mu[j]/=norme;
    }
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul des log-vraisemblances de tous les points*/
/*par rapport a tous les clusters                              */
/*                                                             */
/***************************************************************/

void calculVraisemblancesPointsGaussienneSpherique(int nbIndividus,int nbDimensions,individu_t **individus,
						   int nbClusters,cluster_t *clusters)
{
  /*declaration de variables*/
  int i,j;
  double **points,**pointsRot,**pointsProjetes,**sigma,norme,*angles;
  double logDeterminant,commun;
  /*fin declaration de variables*/

  /*time_t temps_debut,temps_fin,temps_debut2,temps_fin2,temps_coordonnees;
    temps_debut=time(NULL);*/
  
  /*allocation memoire*/
  points=(double **)ckalloc(sizeof(double *)*nbIndividus);
  pointsRot=(double **)ckalloc(sizeof(double *)*nbIndividus);
  pointsProjetes=(double **)ckalloc(sizeof(double *)*nbIndividus);
  for(i=0;i<nbIndividus;i++)
    {
      points[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
      pointsRot[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
      pointsProjetes[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  sigma=(double **)ckalloc(sizeof(double *)*nbDimensions);
  for(i=0;i<nbDimensions;i++)
    {
      sigma[i]=(double *)ckalloc(sizeof(double)*nbDimensions);
    }
  angles=(double *)ckalloc(sizeof(double)*nbDimensions);
  /*fin allocation memoire*/
  
  for(i=0;i<nbIndividus;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  points[i][j]=individus[i]->valeursTraitees[j];
	}
    }
    
  for(i=0;i<nbClusters;i++)
    {
      /*recuperation des angles de mu*/
      cartesiennes2polaires(nbDimensions,clusters[i].mu,angles);
      
      /*temps_debut2=time(NULL);*/
     
      /*calcul des nouvelles coordonnees apres rotation du repere*/
      calculCoordonneesApresRotation(nbDimensions,nbIndividus,points,angles,
				     pointsRot);  
      
      /*calcul des coordonnees des points projetes*/
      calculCoordonneesApresProjection(nbDimensions,nbIndividus,pointsRot,
				       pointsProjetes);
      /*temps_fin2=time(NULL);
	temps_coordonnees+=temps_fin2-temps_debut2;*/

      /*on inverse une copie de la matrice de covariance dans sigma et on
	recupere au passage le log du determinant (regularisee si singuliere)*/
      logDeterminant=inversionMatriceRegularisee(nbDimensions-1,clusters[i].sigma,sigma);
      
      commun=log((double)clusters[i].nbIndividus/(double)nbIndividus)-
	logDeterminant*0.5-(double)(nbDimensions-1)*0.5*log(2.0*M_PI);
      
      /*calcul des vraisemblances*/
      for(j=0;j<nbIndividus;j++)
	{
	  /*individus[j].vrais[i]=-log_determinant*0.5-(double)(nbDimensions-1)*0.5*log(2.0*M_PI);*/
	  norme=calculNorme(nbDimensions-1,pointsProjetes[j],sigma);
	  individus[j]->vrais[i]=-norme*0.5;
	  individus[j]->vrais[i]-=(double)(nbDimensions-1)*log(fabs(1.0+pointsRot[j][0]));
	  individus[j]->vrais[i]+=commun;
	  /*individus[j].vrais[i]+=log((double)clusters[i].nbIndividus/(double)nbIndividus);*/
	}
    }

  /*desallocation memoire*/
  for(i=0;i<nbIndividus;i++)
    {
      ckfree(points[i]);
      ckfree(pointsRot[i]);
      ckfree(pointsProjetes[i]);
    }
  ckfree(points);
  ckfree(pointsRot);
  ckfree(pointsProjetes);
  for(i=0;i<nbDimensions;i++)
    {
      ckfree(sigma[i]);
    }
  ckfree(sigma);
  ckfree(angles);
  /*fin desallocation memoire*/
  
  /*temps_fin=time(NULL);
  printf("vraisemblances : %d\n",temps_fin-temps_debut);
  printf("temps_coordonnees : %d\n",temps_coordonnees);*/
}

/***********************************************/
/*                                             */
/*Procedure de calcul des nouvelles coordonnees*/
/*apres rotation du repere                     */
/*                                             */
/***********************************************/

void calculCoordonneesApresRotation(int nbDimensions,int nbPoints,double **points,
				    double *angles,double **pointsRot)
{
  /*declaration de variables*/
  int i,j;
  double val1,val2,*cosinus,*sinus;
  /*fin declaration de variables*/
  
  /*allocation memoire*/
  cosinus=(double *)ckalloc(sizeof(double)*(nbDimensions-1));
  sinus=(double *)ckalloc(sizeof(double)*(nbDimensions-1));
  /*fin allocation memoire*/

  for(i=0;i<nbDimensions-1;i++)
    {
      cosinus[i]=cos(angles[i]);
      sinus[i]=sin(angles[i]);
    }

  for(i=0;i<nbPoints;i++)
    {
      for(j=0;j<nbDimensions;j++)
	{
	  pointsRot[i][j]=points[i][j];
	}

      for(j=nbDimensions-2;j>=0;j--)
	{
	  val1=pointsRot[i][j];
	  val2=pointsRot[i][j+1];
	
	  pointsRot[i][j]=cosinus[j]*val1+sinus[j]*val2;
	  pointsRot[i][j+1]=-sinus[j]*val1+cosinus[j]*val2;
	}
    }

  /*desallocation memoire*/
  ckfree(cosinus);
  ckfree(sinus);
  /*fin desallocation memoire*/
}

/*********************************************************/
/*                                                       */
/*Procedure de calcul des coordonnees des points projetes*/
/*                                                       */
/*********************************************************/

void calculCoordonneesApresProjection(int nbDimensions,int nbPoints,
				      double **pointsRot,double **pointsProjetes)
{
  /*declaration de variables*/
  int i,j;
  double diviseur;
  /*fin declaration de variables*/

  for(i=0;i<nbPoints;i++)
    {
      if(fabs(1.0+pointsRot[i][0])>ZERO_LIMIT)
	{
	  diviseur=1.0/(1.0+pointsRot[i][0]);
	  for(j=0;j<nbDimensions-1;j++)
	    {
	      pointsProjetes[i][j]=pointsRot[i][j+1]*diviseur;
	    }
	}
      else
	{
	  for(j=0;j<nbDimensions-1;j++)
	    {
	      pointsProjetes[i][j]=9999999.0;
	    }
	}
    }
}


/**********************************************/
/*                                            */
/*Procedure de calcul des coordonnees polaires*/
/*en fonction des coordonnees cartesiennes    */
/*                                            */
/**********************************************/

void cartesiennes2polaires(int nbDimensions,double *cart,double *pol)
{
  /*declaration de variables*/
  int i,j;
  double temp;
  /*fin declaration de variables*/

  pol[0]=acos(cart[0]);
  for(i=1;i<nbDimensions-2;i++)
    {
      pol[i]=cart[i];
      for(j=0;j<i;j++)
	{
	  if(fabs(sin(pol[j]))>ZERO_LIMIT)
	    {
	      pol[i]/=sin(pol[j]);
	    }
	}
      pol[i]=acos(pol[i]);
    }

  pol[nbDimensions-2]=cart[nbDimensions-2];
  temp=cart[nbDimensions-1];
  for(j=0;j<nbDimensions-2;j++)
    {
      if(fabs(sin(pol[j]))>ZERO_LIMIT)
	{
	  pol[nbDimensions-2]/=sin(pol[j]);
	  temp/=sin(pol[j]);
	}
    }

  if(temp>0)
    {
      pol[nbDimensions-2]=acos(pol[nbDimensions-2]);
    }
  else
    {
      pol[nbDimensions-2]=-acos(pol[nbDimensions-2]);
    }
}
