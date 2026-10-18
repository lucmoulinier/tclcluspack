#include "mainClusPack.h"
#include "tools.h"


/******************************************************/
/*                                                    */
/*Procedure de calcul de la distance entre deux points*/
/*                                                    */
/******************************************************/

double calculDistance(int nbDimensions,double *point1,double *point2)
{
  /*declaration des variables*/
  int i;
  double distance,temp;
  /*fin declaration des variables*/

  distance=0;
  for(i=0;i<nbDimensions;i++)
    {
      temp=point1[i]-point2[i];
      distance+=temp*temp;
    }
  distance=sqrt((double)distance);

  return distance;
}

/************************************************************************/
/*                                                                      */
/*Procedure de calcul de la fonction de repartition de la loi de Student*/
/*                                                                      */
/************************************************************************/

double F_Student(double n,double t)
{
  /*declaration des variables*/
  double u;
  /*fin declaration des variables*/
  
  u=(pow(t,2.0/3.0)*(1.0-2.0/(9.0*n))-7.0/9.0)/sqrt(2.0/9.0+pow(t,4.0/3.0)*2.0/(9.0*n));

  return 1.0-F_NormaleCentreeReduite(u);
}

/********************************************************************************/
/*                                                                              */
/*Procedure de calcul de la fonction de repartition de la loi de Fisher-Snedecor*/
/*                                                                              */
/********************************************************************************/

double F_FisherSnedecor(double v1,double v2,double f)
{
  /*declaration des variables*/
  double u;
  /*fin declaration des variables*/
  
  u=(pow(f,1.0/3.0)*(1.0-2.0/(9.0*v2))+2.0/(9.0*v1)-1.0)/
    sqrt(2.0/(9.0*v1)+pow(f,2.0/3.0)*2.0/(9.0*v2));

  return F_NormaleCentreeReduite(u);
}


/*************************************************************************************/
/*                                                                                   */
/*Procedure de calcul de la fonction de repartition de la loi normale_centree_reduite*/
/*                                                                                   */
/*************************************************************************************/

double F_NormaleCentreeReduite(double u)
{
  /*declaration des variables*/
  double b1,b2,b3,b4,b5,t,P;
  /*fin declaration des variables*/

  b1=0.319381530;
  b2=-0.356563782;
  b3=1.781477937;
  b4=-1.821255978;
  b5=1.330274429;
  t=1.0/(1.0+0.2316419*u);

  P=1.0-1.0/sqrt(2.0*M_PI)*exp(-0.5*pow(u,2.0))*(b1*t+b2*pow(t,2.0)+b3*pow(t,3.0)+
						  b4*pow(t,4.0)+b5*pow(t,5.0));

  return P;
}

/************************************************************/
/*                                                          */
/*Procedure de calcul de la fonction de repartition du Chi2,*/
/*(sa probabilite de depasser la valeur x)                  */
/*                                                          */
/************************************************************/

double chiSquare(double x,int k)
{
  /*declaration de variables*/
  int i;
  double resultat;
  /*fin declaration de variables*/

  if(k%2==0)
    {
      /*k pair*/
      resultat=0;
      for(i=0;i<=k/2-1;i++)
	{
	  resultat+=pow(x/2.0,(double)i)/(double)calculFactorielle(i);
	}
      resultat*=exp(-x/2.0);
    }
  else
    {
      /*k impaire*/
      resultat=0;
      for(i=1;i<=(k-1)/2;i++)
	{
	  resultat+=pow(x,(double)i-0.5)/(double)calculFactorielleImpaire(i);
	}
      resultat=2.0*(1.0-F_NormaleCentreeReduite(sqrt(x)))+
	2.0/sqrt(2.0*M_PI)*exp(-x/2.0)*resultat;
    }

  return resultat;
}


/***************************************************************************/
/*                                                                         */
/*Procedure de regression lineaire sur les mesures (renvoie le coefficient)*/
/*                                                                         */
/***************************************************************************/

double regressionLineaire(int n,double **mesures)
{
  /*declaration des variables*/
  int i;
  double moyennex,moyenney,covariance,variancex;
  /*fin declaration des variables*/

  /*calcul des moyennes*/
  moyennex=0;
  moyenney=0;

  for(i=0;i<n;i++)
    {
      moyennex+=mesures[i][0];
      moyenney+=mesures[i][1];
    }
  
  moyennex/=(double)n;
  moyenney/=(double)n;

  /*calcul de la variance de x*/
  variancex=0;

  for(i=0;i<n;i++)
    {
      variancex+=pow((double)(moyennex-mesures[i][0]),2.0);
    }
  variancex/=(double)n;

  /*calcul de la covariance entre x et y*/
  covariance=0;
  
  for(i=0;i<n;i++)
    {
      covariance+=(mesures[i][0]-moyennex)*(mesures[i][1]-moyenney);
    }
  covariance/=(double)n;

  return covariance/variancex;
}

/*****************************************************/
/*                                                   */
/*renvoie la valeur absolue de a                     */
/*                                                   */
/*****************************************************/
double valeurAbsolue(double a)
{
  if(a>0)
    {
      return a;
    }
  else
    {
      return -a;
    }
}

 
/*****************************************************/
/*                                                   */
/*renvoie le minimum de a et de b                    */
/*                                                   */
/*****************************************************/
int min(int a,int b)
{
  if(a<b)
    {
      return a;
    }
  else
    {
      return b;
    }
}

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapide(double *valeurs,int gauche,int droite)
{ 
  /*declaration des variables*/
  int elementSuivant;
  int indiceSeparateur;
  /*fin de declaration des variables*/
  
  if(gauche>=droite)
    {
      return;
    }
  
  indiceSeparateur=gauche;
  elementSuivant=gauche+1;
  
  while(elementSuivant<=droite)
    {
      if(valeurs[elementSuivant]>=valeurs[indiceSeparateur])
	{
	  echangerValeursPourTri(valeurs,indiceSeparateur+1,elementSuivant);
	  echangerValeursPourTri(valeurs,indiceSeparateur,indiceSeparateur+1);
	  indiceSeparateur++;
	}
      elementSuivant++;
    }
  
  triRapide(valeurs,gauche,indiceSeparateur-1);
  triRapide(valeurs,indiceSeparateur+1,droite);
}

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echangerValeursPourTri(double *valeurs,int element1,int element2)
{
  /*declaration des variables*/
  double tempVal;
  /*fin de declaration des variables*/

  tempVal=valeurs[element1];
  valeurs[element1]=valeurs[element2];
  valeurs[element2]=tempVal;
}

/****************************************/
/*                                      */
/*Procedure de tri par ordre decroissant*/
/*                                      */
/****************************************/

void triRapideValeursPositions(double *valeurs,int *positions,int gauche,int droite)
{ 
  /*declaration des variables*/
  int elementSuivant;
  int indiceSeparateur;
  /*fin de declaration des variables*/
  
  if(gauche>=droite)
    {
      return;
    }
  
  indiceSeparateur=gauche;
  elementSuivant=gauche+1;
  
  while(elementSuivant<=droite)
    {
      if(valeurs[elementSuivant]>=valeurs[indiceSeparateur])
	{
	  echangerValeursPositionsPourTri(valeurs,positions,indiceSeparateur+1,
					  elementSuivant);
	  echangerValeursPositionsPourTri(valeurs,positions,indiceSeparateur,
					  indiceSeparateur+1);
	  indiceSeparateur++;
	}
      elementSuivant++;
    }
  
  triRapideValeursPositions(valeurs,positions,gauche,indiceSeparateur-1);
  triRapideValeursPositions(valeurs,positions,indiceSeparateur+1,droite);
}

/******************************/
/*                            */
/*sous-procedure de tri_rapide*/ 
/*                            */
/******************************/

void echangerValeursPositionsPourTri(double *valeurs,int *positions,int element1,int element2)
{
  /*declaration des variables*/
  int tempPosition;
  double tempVal;
  /*fin de declaration des variables*/

  tempVal=valeurs[element1];
  valeurs[element1]=valeurs[element2];
  valeurs[element2]=tempVal;

  tempPosition=positions[element1];
  positions[element1]=positions[element2];
  positions[element2]=tempPosition;
}

/****************************************************************/
/*                                                              */
/*Procedure de test d'ajustement avec la distribution gaussienne*/
/*par Smirnov-Kolmogorov                                        */
/*                                                              */
/****************************************************************/

int testAjustementDistributionGaussienne(int nbValeurs,double *valeurs,double moyenne,
					    double ecartType)
{
  /*declaration de variables*/
  int i;
  double Dn,differenceCourante,seuil;
  double valeursLimiteRisque_01[100]=
  {0.995, 0.92929, 0.829, 0.73424, 0.66853, 0.61661, 0.57581, 0.54179, 0.51332, 0.48893,
   0.4677, 0.44905, 0.43247, 0.41762, 0.40420, 0.39201, 0.38086, 0.37062, 0.36117, 0.35241,
   0.34427, 0.33666, 0.32954, 0.32286, 0.31657, 0.31064, 0.30502, 0.29971, 0.29466, 0.28987,
   0.2853, 0.28094, 0.27677, 0.27279, 0.26897, 0.26532, 0.26180, 0.25843, 0.25518, 0.25205,
   0.24904, 0.24613, 0.24332, 0.2406, 0.23798, 0.23544, 0.23298, 0.23059, 0.22828, 0.22604,
   0.22386, 0.22174, 0.21968, 0.21768, 0.21574, 0.21384, 0.21199, 0.21019, 0.20844, 0.20673,
   0.20506, 0.20343, 0.20184, 0.20029, 0.19877, 0.19729, 0.19584, 0.19442, 0.19303, 0.19167,
   0.19034, 0.18903, 0.18776, 0.1865, 0.18528, 0.18408, 0.18290, 0.18174, 0.18060, 0.17949,
   0.1784, 0.17732, 0.17627, 0.17523, 0.17421, 0.17321, 0.17223, 0.17126, 0.17031, 0.16938,
   0.16846, 0.16755, 0.16666, 0.16579, 0.16493, 0.16408, 0.16324, 0.16242, 0.16161, 0.16081};
  /*fin declaration de variables*/
  printf("debut\n");
  Dn=0;
  for(i=0;i<nbValeurs;i++)
    {
      if(valeurs[i]>moyenne)
	{
	  break;
	}
      differenceCourante=fabs(F_NormaleCentreeReduite((valeurs[i]-moyenne)/ecartType)-
			       ((double)nbValeurs-(double)i)/(double)nbValeurs);
      if(differenceCourante>Dn)
	{
	  Dn=differenceCourante;
	}
    }

  if(nbValeurs>100)
    {
      seuil=1.629/sqrt((double)nbValeurs);
    }
  else
    {
      seuil=valeursLimiteRisque_01[nbValeurs-1];
    }
  printf("seuil: %f\n",seuil);
  printf("Dn: %f\n",Dn);

printf("fin\n");
  if(Dn>seuil)
    {
      /*pas d'ajustement*/
      return NON;
    }
  else
    {
      /*il y a ajustement*/
      return OUI;
    }
}

/*****************************************************/
/*                                                   */
/*Procedure de calcul de la factorielle d'un entier n*/
/*                                                   */
/*****************************************************/

int calculFactorielle(int n)
{
  /*declaration de variables*/
  int i,resultat;
  /*fin declaration de variables*/

  resultat=1;
  for(i=2;i<=n;i++)
    {
      resultat*=i;
    }
  return resultat;
}

/***************************************************************/
/*                                                             */
/*Procedure de calcul de la factorielle "impaire" d'un entier n*/
/*                                                             */
/***************************************************************/

int calculFactorielleImpaire(int n)
{
  /*declaration de variables*/
  int i,resultat;
  /*fin declaration de variables*/

  resultat=1;
  for(i=2;i<=n;i++)
    {
      resultat*=2*i-1;
    }
  return resultat;
}
