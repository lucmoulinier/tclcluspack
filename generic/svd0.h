#define SIGN(a,b) ((b) >= 0.0 ? ABS(a) : -ABS(a))

double ABS(double a);

double DMAX(double a,double b);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void sswap(double **v, int n, int i, int id);

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
/*	by : Doulaye Demebele                                                    */
/*---------------------------------------------------------------------------*/
void ce2svd(double **A, int m, int n, double **U, double *W, 
	    double **V, int *ierr);
