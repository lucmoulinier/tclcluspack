#include "mainClusPack.h"
#include "svd0.h"

double ABS(double a)
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


double DMAX(double a,double b)
{
  if(a>b)
    {
      return a;
    }
  else
    {
      return b;
    }
}



/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void sswap(double **v, int n, int i, int id)
{	int ii, m, mp1;
  double tmp, *sx, *sy;

  /*allocation memoire*/
  sx = (double *)ckalloc(sizeof(double)*(n+1)); 
  sy = (double *)ckalloc(sizeof(double)*(n+1));
  /*fin allocation memoire*/

  for (ii=1; ii<=n; ii++)	sx[ii] = v[ii][i];
  for (ii=1; ii<=n; ii++)	sy[ii] = v[ii][id];
  m = n%3;
  if (m==0) goto c40;
  for (ii=1; ii<=m; ii++)
    {	tmp = sx[ii];	sx[ii] = sy[ii]; sy[ii] = tmp;
    }
  if (n<3) goto c50;
 c40:;
  mp1 = m+1;
  for (ii=mp1; ii<=n; ii += 3)
    {	tmp = sx[ii];	sx[ii]   = sy[ii];   sy[ii]   = tmp;
      tmp = sx[ii+1];	sx[ii+1] = sy[ii+1]; sy[ii+1] = tmp;
      tmp = sx[ii+2];	sx[ii+2] = sy[ii+2]; sy[ii+2] = tmp;
    }
 c50:;
	
  for (ii=1; ii<=n; ii++)	v[ii][i] = sx[ii];
  for (ii=1; ii<=n; ii++)	v[ii][id] = sy[ii];

  /*desallocation memoire*/
  ckfree(sx);
  ckfree(sy);
  /*fin desallocation memoire*/
}

/*---------------------------------------------------------------------------*/
/*	Compute the Singular Value Decomposition of matrix A(m,n), where  m >= n */
/*	                    A = U(m,n)W(n,n)V^T(n,n)                             */
/*	if m is less than n, ierr is set to -1 and no calculations are done      */
/*                                                                           */
/*	Translated from the FORTRAN function which appeared in :                 */
/*	ACM Transaction on Mathematics Software, vol.8, n.1, mar.1982, pp.84-    */
/*	by Tony Chan.                                                            */
/*	The original function is an ALGOL one by Golub & Reinsch which appeared  */
/*	in Num. Math. vol.14, pp.403-420, (1970) ; see also                      */
/*	Handbook for Automatic Comp. vol.II Linear Algebra, pp.134-151 (1971)    */
/*                                                                           */
/*	Date : 28 October 1998                                                   */
/*	by : Doulaye Dembele                                                    */
/*---------------------------------------------------------------------------*/

void ce2svd(double **A, int m, int n, double **U, double *W, double **V, int *ierr)
{
  int i, i1, ii, its, j, k, k1, kk, l, l1, ll, mn, nm1, ierrp1, id, ip1;
  double c, f, g, h, x, s, y, z, scale, *rv1, eps, **u, **v, *w;
  
  /*allocation memoire*/
  u=(double **)ckalloc(sizeof(double *)*(m+1));
  for(i=0;i<m+1;i++)
    {
      u[i]=(double *)ckalloc(sizeof(double)*(n+1));
    }
  v=(double **)ckalloc(sizeof(double *)*(n+1));
  for(i=0;i<n+1;i++)
    {
      v[i]=(double *)ckalloc(sizeof(double)*(n+1));
      for(j=0;j<n+1;j++)
	{
	  v[i][j]=0.0000001;
	}
    }
  w=(double *)ckalloc(sizeof(double)*(n+1));
  rv1=(double *)ckalloc(sizeof(double)*(m+1));
  /*fin allocation memoire*/
  
  for (i=0; i<m; i++) for (j=0; j<n; j++) 
			{
			  u[i+1][j+1]=A[i][j];
			}
  
  *ierr = 0;
  if (m<n) {*ierr = -1; goto c670; }

  /* Householder reduction to bidiagonal form */
  g = scale = x = 0.0;
  for (i=1; i<=n; i++)
    {
      l = i+1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
		
      /* Compute left transformation that zero the subdiagonal elements
	 of the i-th column */
      for (k=i; k<=m; k++)	scale += ABS(u[k][i]);
      if (scale)
	{	for (k=i; k<=m; k++)	{ u[k][i] /= scale; s += u[k][i]*u[k][i]; }
	  f = u[i][i];
	  g = -SIGN(sqrt(s),f);
	  h = f*g-s;
	  u[i][i] = f-g;
	  if (i<n)
	    {	/* Apply left transformation to remaining columns of A */
	      for (j=l; j<=n; j++)
		{	s = 0.0;
		  for (k=i; k<=m; k++)	s += u[k][i]*u[k][j];
		  f = s/h;
		  for (k=i; k<=m; k++)	u[k][j] += f*u[k][i];
		}
	    }
	  /* Compute right transformations */
	  for (k=i; k<=m; k++)	u[k][i] *= scale;
	}
      w[i] = scale*g;
      g = s = scale = 0.0;
      if ((i<=m) || (i!=n))
	{	for (k=l; k<=n; k++)	scale += ABS(u[i][k]);
	  if (scale)
	    {	for (k=l; k<=n; k++)	
		{	u[i][k] /= scale; 
		  s += u[i][k]*u[i][k]; 
		}
	      f = u[i][l];
	      g = -SIGN(sqrt(s),f);
	      h = f*g-s;
	      u[i][l] = f-g;
	      for (k=l; k<=n; k++)	rv1[k] = u[i][k]/h;
	      if (i!=m)
		{	for (j=l; j<=m; j++)
		    {	s = 0.0;
		      for (k=l; k<=n; k++)	s += u[j][k]*u[i][k];
		      for (k=l; k<=n; k++)	u[j][k] += s*rv1[k];
		    }
		}
	      for (k=l; k<=n; k++)	u[i][k] *= scale;
	    }
	}
      x = DMAX(x,ABS(w[i])+ABS(rv1[i]));
    }

  /* Accumulation of the right-hand transformations */
  for (ii=1; ii<=n; ii++)
    {	i = n+1-ii;
      if (i<n)
	{	if (g)
	    {	/* Double division avoids possible underflow */
	      for (j=l; j<=n; j++)	v[j][i] = (u[i][j]/u[i][l])/g;
	      for (j=l; j<=n; j++)
		{	s = 0.0;
		  for (k=l; k<=n; k++)	s += u[i][k]*v[k][j];
		  for (k=l; k<=n; k++)	v[k][j] += s*v[k][i];
		}
	    }
	  for (j=l; j<=n; j++)	{ v[i][j] = v[j][i] = 0.0; }
	}
      v[i][i] = 1.0;	g = rv1[i]; l = i;
    }
	
  /* Accumulation of the left-hand transformations */
  /* For i=min(m,n)-1 step -1 until 0 do */
  mn = n;
  if (m<n) mn = m;
  for (ii=1; ii<=mn; ii++)
    {	i = mn+1-ii;	l = i+1;	g = w[i];
      if (i!=n)
	for (j=l; j<=n; j++)	u[i][j] = 0.0;
      if (g)
	{	if (i!=mn)
	    {	for (j=l; j<=n; j++)
		{	s = 0.0;
		  for (k=l; k<=m; k++)	s += u[k][i]*u[k][j];
			
		  /* Double division avoids possible underflow */
		  f = (s/u[i][i])/g;
		  for (k=i; k<=m; k++)	u[k][j] += f*u[k][i];
		}
	    }
	  for (j=i; j<=m; j++)	u[j][i] /= g;
	  goto c450;
	}
      for (j=i; j<=m; j++)	u[j][i] = 0.0;
    c450:;
      u[i][i] += 1.0;
    }
	
  /* Diagonalization of the bidagonal form */

  /*Modif Nico	eps = 2.22e-16*x;*/
  eps=0.00000001;
  /* For k=n step -1 until 1 do */
  for (kk=1; kk<=n; kk++)
    {	k1 = n-kk;	k = k1+1; its = 0;
		
      /* Test for splitting : for k=k step -1 until 1 do */
    c480:;
      for (ll=1; ll<=k; ll++)
	{	l1 = k-ll;  l = l1+1;
	  if (ABS(rv1[l])<=eps) goto c550;
			
	  /* rv1(1) is always zero, so there is no exit through 
	     the bottom of the loop*/
	  if (ABS(w[l1])<=eps) goto c500;
	}

      /* Cancellation of rv1(l) if l greater than 1*/
    c500:;
      c = 0.0;	s = 1.0;
      for (i=l; i<=k; i++)
	{	f = s*rv1[i];
	  rv1[i] = c*rv1[i];
	  if (ABS(f)<=eps) goto c550;
	  g = w[i];
	  h = sqrt(f*f+g*g);
	  w[i] = h;
	  c = g/h;
	  s = -f/h;
			
	  /*------- fmodification of matrix u --------*/
	  for (j=1; j<=m; j++)
	    {	y = u[j][l1];
	      z = u[j][i];
	      u[j][l1] = y*c+z*s;
	      u[j][i] = -y*s+z*c;
	    }
	}

      /* Test for convergence */
    c550:;
      z = w[k];
      if (l!=k)
	{	/* Shift from bottom 2 by 2 minor */
	  if (its==30) goto c660;
	  its++;
	  x = w[l];	y = w[k1]; g = rv1[k1]; h = rv1[k]; 
	  f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g = sqrt(f*f+1.0);
	  f = ((x-z)*(x+z)+h*(y/(f+SIGN(g,f))-h))/x;
				
	  /* Next QR transformation*/
	  c = s = 1.0;		
	  for (i1=l; i1<=k1; i1++)
	    {	i = i1+1;
	      g = rv1[i];	y = w[i]; h = s*g; g = c*g; z = sqrt(f*f+h*h); 
	      rv1[i1]= z; c = f/z; s = h/z; f = x*c+g*s; g = -x*s+g*c; 
	      h = y*s; y = y*c;

	      /*---- modification of matrix v ------------*/
	      for (j=1; j<=n; j++)
		{	x = v[j][i1]; z = v[j][i]; 
		  v[j][i1] = x*c+z*s; v[j][i] = -x*s+z*c;
		}
	      z = sqrt(f*f+h*h);
	      w[i1] = z;
					
	      /* Rotation can be arbitrary if z is zero */
	      if (z) {	c = f/z; s = h/z; }
	      f = c*g+s*y; x = -s*g+c*y;

	      /*---- modification of matrix u ------------*/
	      for (j=1; j<=m; j++)
		{	y = u[j][i1]; z = u[j][i];
		  u[j][i1] = y*c+z*s; u[j][i] = -y*s+z*c;
		}
	    }
	  rv1[l] = 0.0;
	  rv1[k] = f;
	  w[k] = x;
	  goto c480;
	}
      /* Convergence */
      if (z<0.0)
	{	/* w(k) is made non-negative*/
	  w[k] = -z;

	  /*---- modification of matrix v ------------*/
	  for (j=1; j<=n; j++)	v[j][k] = -v[j][k];
	}
    }
  goto c670;
  /* Set error - 
     no convergence to a singular value after 30 iterations */
 c660:;
  *ierr = k;
 c670:;

  /*desallocation memoire*/
  ckfree(rv1);
  /*fin desallocation memoire*/

  ierrp1 = *ierr+1;
  if ((*ierr>=0) && (n>1) && (ierrp1!=n))
    {	/*	Sort singular values and exchange columns of U and V accordingly */
      /*	Selection sort minimizes swapping of U and V */
      nm1 = n-1;
      for (i=ierrp1; i<=nm1; i++)
	{	/*	Find index of maximum singular value */
	  id = i;	ip1 = i+1;
	  for (j=ip1; j<=n; j++)	if (w[j]>w[id])	id = j;
	  if (id!=i)
	    {	/*	Swap singular values and vectors */
	      x = w[i]; w[i] = w[id];	w[id] = x;
	      /*	Swap v(1,i) and v(1,id) */
	      sswap(v,n,i,id);
	      /*	Swap u(1,i) and u(1,id) */
	      sswap(u,m,i,id);
	    }
	}
    }
  /*	Now come back to matrices U, V and vector W */
  for (i=0; i<m; i++)	for (j=0; j<n; j++)	U[i][j] = u[i+1][j+1];
  for (i=0; i<n; i++)	for (j=0; j<n; j++)	V[i][j] = v[i+1][j+1];
  for (j=0; j<n; j++)	W[j] = w[j+1];

  /*desallocation memoire*/
  ckfree(w);
  for(i=0;i<n+1;i++)
    {
      ckfree(v[i]);
    }
  ckfree(v);
  for(i=0;i<m+1;i++)
    {
      ckfree(u[i]);
    }
  ckfree(u);
  /*fin desallocation memoire*/

}/*-------------------------- end of function ce2svd() ----------------------*/













