/*
 * File name: dslab.c
 * Date:      2016/06/19 18:10
 * Author:    Jiri Brozovsky

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; version 2 of the
   License.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   in a file called COPYING along with this program; if not, write to
   the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
   02139, USA.

 */

#include <stdio.h>
#ifdef LARGE
#include <stdlib.h>
#ifdef UI
#include <string.h>
#endif
#endif
#include <math.h>

int    n_nodes = 0 ;
int    n_elems = 0 ;
int    n_disps = 0 ;
int    n_nfors = 0 ;

double *x_i = NULL ;
double *y_i = NULL ;

double *E = NULL ;
double *nu = NULL ;
double *tl = NULL ;
int   *n1 = NULL ;
int   *n2 = NULL ;
int   *n3 = NULL ;
int   *n4 = NULL ;
int   *e_g = NULL ; /* element group */

int   *d_n = NULL ; /* node */
int   *d_d = NULL ; /* direction 1=w 2=rotx 3=roty */
double *d_v = NULL ; /* size */
int   *d_g = NULL ; /* load group */

int   *f_n = NULL ; /* node */
int   *f_d = NULL ; /* direction 1=fx 2=fy 3=m */
double *f_v = NULL ; /* size */
int   *f_g = NULL ; /* load group */

/* solution variables: */
double  ke[12][12] ;
double  B[3][12] ;
double  SB[12][3] ;
double  SBD[12][3] ;
double  SBDB[12][12] ;
double  SBDBS[12][12] ;
double  D[3][3] ;
double  fe[12];
double  ue[12];
double  **S;

int    K_len    = 0 ;
int    K_size   = 0 ;
int   *K_sizes  = NULL ; /* lenghts of K's rows */
int   *K_from   = NULL ;
int   *K_cols   = NULL ;
double *K_val    = NULL ;
double *F_val    = NULL ;
double *u_val    = NULL ;

double *M    = NULL ;
double *r    = NULL ;
double *z    = NULL ;
double *p    = NULL ;
double *q    = NULL ;

/** Reading of data from file */
int read_data(fw)
FILE *fw ;
{
  int i;

  fprintf(stderr,"Number of nodes:\n");
  fscanf(fw,"%d", &n_nodes);
  if (n_nodes <= 0)
  {
    fprintf(stderr,"No nodes!\n");
    return(-1);
  }

  if ((S = (double **)malloc(12*sizeof(double *))) != NULL)
  {
    for (i=0; i<12; i++) S[i] = (double *)malloc(12*sizeof(double)) ;
  } else { exit(-1); }

  /* allocate data for nodes */
  if ((x_i=(double *)malloc(n_nodes*sizeof(double))) == NULL)
  {
    fprintf(stderr,"Can't allocate X coordinates!\n");
    return(-2);
  }

  if ((y_i=(double *)malloc(n_nodes*sizeof(double))) == NULL)
  {
    free(x_i);
    fprintf(stderr,"Can't allocate Y coordinates!\n");
    return(-2);
  }

  /* read coordinates of nodes */
#ifdef UI
  fprintf(stderr,"Coordinates of nodes (x y):\n");
#endif
  for (i=0; i<n_nodes; i++)
  {
    fscanf(fw,"%lf %lf", &x_i[i], &y_i[i]) ;
  }
#ifdef UI
  fprintf(stderr,"  Have %d coordinates.\n",n_nodes);
#else
  fprintf(stderr,"  Nodes:           %d\n",n_nodes);
#endif


  /* number of elements */
#ifdef UI
  fprintf(stderr,"Number of elements:\n");
#endif
  fscanf(fw,"%d", &n_elems);
  if (n_elems <= 0)
  {
    free(x_i); free(y_i); 
    fprintf(stderr,"No elements!\n");
    return(-1);
  }

  if ((n1=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); 
    fprintf(stderr,"Can't allocate node positions!\n");
    return(-2);
  }
  if ((n2=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1);
    fprintf(stderr,"Can't allocate node positions!\n");
    return(-2);
  }
  if ((n3=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);
    fprintf(stderr,"Can't allocate node positions!\n");
    return(-2);
  }
  if ((n4=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3);
    fprintf(stderr,"Can't allocate node positions!\n");
    return(-2);
  }

  if ((E=(double *)malloc(n_elems*sizeof(double))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);
    fprintf(stderr,"Can't allocate stifnesses!\n");
    return(-2);
  }
  if ((nu=(double *)malloc(n_elems*sizeof(double))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);
    fprintf(stderr,"Can't allocate areas!\n");
    return(-2);
  }

  if ((tl=(double *)malloc(n_elems*sizeof(double))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);
    fprintf(stderr,"Can't allocate areas!\n");
    return(-2);
  }

  if ((e_g=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);
    fprintf(stderr,"Can't allocate types!\n");
    return(-2);
  }


  /* read element data */
#ifdef UI
  fprintf(stderr,"Element data (node1 node2 node3 node4 E nu thick grp):\n");
#endif
  for (i=0; i<n_elems; i++)
  {
    fscanf(fw,"%d %d %d %d %lf %lf %lf %d",
      &n1[i],&n2[i],&n3[i],&n4[i],&E[i],&nu[i],&tl[i],&e_g[i]) ;
  }
#ifdef UI
  fprintf(stderr,"  Have %d elements.\n",n_elems);
#else
  fprintf(stderr,"  Elements:        %d\n", n_elems);
#endif

  /* supports */
#ifdef UI
  fprintf(stderr,"Number of supports:\n");
#endif
  fscanf(fw,"%d", &n_disps);
  if (n_disps <= 0)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    fprintf(stderr,"No supports!\n");
    return(-1);
  }

  if ((d_n=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    fprintf(stderr,"Can't allocate nodes for supports!\n");
    return(-2);
  }
  if ((d_d=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);
    fprintf(stderr,"Can't allocate types of supports!\n");
    return(-2);
  }
  if ((d_v=(double *)malloc(n_disps*sizeof(double))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);free(d_d);
    fprintf(stderr,"Can't allocate sizes of supports!\n");
    return(-2);
  }
  if ((d_g=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);free(d_d);free(d_v);
    fprintf(stderr,"Can't allocate load groups!\n");
    return(-2);
  }

  
  /* read supports data */
#ifdef UI
  fprintf(stderr,"Supports data (node direction size grp):\n");
#endif
  for (i=0; i<n_disps; i++)
  {
    fscanf(fw,"%d %d %lf %d",&d_n[i], &d_d[i], &d_v[i], &d_g[i]) ;
  }
#ifdef UI
  fprintf(stderr,"  Have %d supports.\n",n_disps);
#else
  fprintf(stderr,"  Supports:        %d\n",n_disps);
#endif


  /* forces in nodes */
#ifdef UI
  fprintf(stderr,"Number of forces in nodes:\n");
#endif
  fscanf(fw,"%d", &n_nfors);
  if (n_nfors <= 0)
  {
    fprintf(stderr,"No forces in nodes!\n");
  }
  else
  {
  if ((f_n=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    fprintf(stderr,"Can't allocate nodes for forces!\n");
    return(-2);
  }
  if ((f_d=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    free(f_n);
    fprintf(stderr,"Can't allocate types of forces!\n");
    return(-2);
  }
  if ((f_v=(double *)malloc(n_nfors*sizeof(double))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    free(f_n);free(f_d);
    fprintf(stderr,"Can't allocate sizes of forces!\n");
    return(-2);
  }
  if ((f_g=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    free(f_n);free(f_d);free(f_v);
    fprintf(stderr,"Can't allocate load groups!\n");
    return(-2);
  }

  
  /* read forces data */
#ifdef UI
  fprintf(stderr,"Forces in nodes data (node direction size grp):\n");
#endif
  for (i=0; i<n_nfors; i++)
  {
    fscanf(fw,"%d %d %lf %d",&f_n[i], &f_d[i], &f_v[i], &f_g[i]) ;
  }
#ifdef UI
  fprintf(stderr,"  Have %d forces in nodes.\n",n_nfors);
#else
  fprintf(stderr,"  Forces in nodes: %d\n",n_nfors);
#endif
  } /* end of n_nfors > 0 */

  fprintf(stderr,"End of input.\n");

  return(0);
}

/** Write input data to file */
#ifdef UI
int write_idata(fw)
FILE *fw ;
{
  int i;

  fprintf(fw,"%d\n", n_nodes);
  for (i=0; i<n_nodes; i++)
  { fprintf(fw,"%e %e\n", x_i[i], y_i[i]) ; }

  fprintf(fw,"%d\n", n_elems);
  for (i=0; i<n_elems; i++)
  {
    fprintf(fw,"%d %d %d %d %e %e %e %d\n",
      n1[i],n2[i],n3[i],n4[i],E[i],nu[i],tl[i],e_g[i]) ;
  }

  fprintf(fw,"%d\n", n_disps);
  for (i=0; i<n_disps; i++)
  { fprintf(fw,"%d %d %e %d\n",d_n[i], d_d[i], d_v[i], d_g[i]) ; }

  fprintf(fw,"%d\n", n_nfors);
  for (i=0; i<n_nfors; i++)
  { fprintf(fw,"%d %d %e %d\n",f_n[i], f_d[i], f_v[i], f_g[i]) ; }

  return(0);
}
#endif

/** free K,F */
void free_sol_data()
{
  if (F_val != NULL )free(F_val);
  if (u_val != NULL )free(u_val);
  if (K_sizes != NULL )free(K_sizes);
  if (K_from != NULL )free(K_from);
  if (K_cols != NULL )free(K_cols);
  if (K_val != NULL )free(K_val);

  if (M != NULL )free(M);
  if (r != NULL )free(r);
  if (z != NULL )free(z);
  if (p != NULL )free(p);
  if (q != NULL )free(q);
  
}


/* Allocates space for linear system */
int alloc_kf()
{
  int i,j, sum;
  int k_size = 0 ;

  k_size = 3*n_nodes ;

  if ((K_sizes = (int *)malloc(k_size*sizeof(int)))   == NULL) { goto memFree;}
  if ((K_from  = (int *)malloc(k_size*sizeof(int)))   == NULL) {goto memFree;} 
  if ((F_val = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 
  if ((u_val = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 

  if ((M   = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 
  if ((r   = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 
  if ((z   = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 
  if ((p   = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 
  if ((q   = (double *)malloc(k_size*sizeof(double))) == NULL) {goto memFree;} 

  K_len = 0 ;

  for (i=0; i<k_size; i++) 
  { 
    K_sizes[i] = 0 ; 
    K_from[i]  = 0 ; 
    F_val[i]   = 0.0 ; 
    u_val[i]   = 0.0 ; 
    M[i]       = 0.0 ; 
    r[i]       = 0.0 ; 
    z[i]       = 0.0 ; 
    p[i]       = 0.0 ; 
    q[i]       = 0.0 ; 
  } 

  for (i=0; i<n_nodes; i++)
  {
    for (j=0; j<n_elems; j++)
    {
      if ((n1[j]-1) == i) 
      {
        K_sizes[3*i+0]+=12; 
        K_sizes[3*i+1]+=12; 
        K_sizes[3*i+2]+=12; 
      }
      if ((n2[j]-1) == i) 
      {
        K_sizes[3*i+0]+=12; 
        K_sizes[3*i+1]+=12; 
        K_sizes[3*i+2]+=12; 
      }
      if ((n3[j]-1) == i) 
      {
        K_sizes[3*i+0]+=12; 
        K_sizes[3*i+1]+=12; 
        K_sizes[3*i+2]+=12; 
      }
      if ((n4[j]-1) == i) 
      {
        K_sizes[3*i+0]+=12; 
        K_sizes[3*i+1]+=12; 
        K_sizes[3*i+2]+=12; 
      }
    }
  }

  sum = 0 ;

  for (i=0; i<k_size; i++)
  {
    K_from[i] = sum ;
    sum += K_sizes[i] ;
  }

  K_len = sum ;

  if ((K_cols = (int *)malloc(K_len*sizeof(int))) == NULL) { goto memFree ; } 
  if ((K_val = (double *)malloc(K_len*sizeof(double))) == NULL) { goto memFree ; } 

  for (i=0; i<K_len; i++)
  {
    K_cols[i] = -1 ;
    K_val[i]  = 0.0 ;
  }

  K_size = k_size ;
  
  return(0);
memFree:
  fprintf(stderr,"Not enough memory!");
  free_sol_data();
  return(-1);
}

double norm_K()
{
  int i,j ;
  double MaxNorm = 0.0 ;
  double Norm    = 0.0 ;
  
  for (i=0; i<(3*n_nodes); i++)
  {
     Norm = 0.0;
    for (j=K_from[i]; j<K_from[i]+K_sizes[i]; j++)
    {
      if (K_cols[j] < 0) {break;}
      Norm += (K_val[j]*K_val[j]);
    }
    Norm = sqrt(Norm);
    if (Norm > MaxNorm) {MaxNorm = Norm;}
  }
  return(Norm);
}

double vec_norm(a, len)
double *a;
int len;
{
  int i ;
  double Norm    = 0.0 ;
  
  for (i=0; i<len; i++) { Norm += (a[i]*a[i]); }
  return(sqrt(Norm));
}

int solve_eqs()
{
  double ro, alpha, beta;
  double roro = 0.0 ;
  double normRes, normX, normA, normB;
  double mval ;
  int   converged = 0;
  int   n = 0;
  int   i,j,k;

  n = K_size;

  normA = norm_K();
  normB = vec_norm(F_val, n);

  if (normB <= 0.0) /* no loads - nothing to do */
  {
    fprintf(stderr,"No load found!\n");
    return(0);
  }

  /* Jacobi preconditioner: */
  for (i=0; i<n; i++) 
  { 
    M[i] = 0.0 ;
    for (j=K_from[i]; j<K_from[i]+K_sizes[i]; j++)
    {
      if (K_cols[j] == (i))
      {
        M[i] = K_val[j] ;
        break ;
      }
    }

    if (fabs(M[i]) < 1e-5) 
    { 
      fprintf(stderr,"zero value at [%d,%d]: %e\n",i+1,i+1, M[i]);
      return( -1 ); 
    }
  }

  /* r = b - A*x  */
  for (i=0; i<n; i++)
  {
    mval = 0.0 ;

    for (j=0; j<K_sizes[i]; j++)
    {
      if  (K_cols[K_from[i]+j] < 0) {break;}
      mval += K_val[K_from[i]+j] * u_val[K_cols[K_from[i]+j]];
    }
    r[i] = mval ;
  }
  for (i=0; i<n; i++) { r[i] = F_val[i] - r[i] ; }

  /* main loop */
  for (i=1; i<=n; i++) 
  { 
    fprintf(stderr,"  CG iteration: %d/%d\n",i,n);
    for (j=0; j<n; j++) { z[j] = (r[j] / M[j]) ; }

    ro = 0.0 ;
    for (j=0; j<n; j++) {ro += r[j]*z[j];}

    if (i == 1)
    {
      for (j=0; j<n; j++) { p[j] = z[j]; }
    }
    else
    {
      beta = ro / roro ;
      for (j=0; j<n; j++) { p[j] = (z[j] + (beta*p[j])) ; }
    }

    for (k=0; k<n; k++) /* q = K*p */
    {
      mval = 0.0 ;

      for (j=0; j<K_sizes[k]; j++)
      {
        if  (K_cols[K_from[k]+j] < 0) {break;}
        mval += K_val[K_from[k]+j] * p[K_cols[K_from[k]+j]];
      }
      q[k] = mval ;
    }

    mval = 0.0 ;
    for (j=0; j<n; j++) {mval += p[j]*q[j];}
    alpha = ro / mval ;

    for (j=0; j<n; j++) 
    { 
      u_val[j] = u_val[j] + (alpha * p[j])  ; 
      r[j] = r[j] - (alpha * q[j])  ; 
    } 

    /* Convergence testing */
    normRes = vec_norm(r, n);
    normX   = vec_norm(u_val, n);

    if (normRes  <= ((1e-3)*((normA*normX) + normB)) ) 
    {
      converged = 1;
      break;
    }

    roro = ro;
  
  } /* end of main loop */

  if (converged == 1) { return(0); }
  else                
  {
    fprintf(stderr,"Unconverged solution!\n");
    return(-1); 
  }
}



/* puts data to the right place in K */
void md_K_add(row, col, val)
int row;
int col;
double val;
{
  int i ;

  for (i=K_from[row-1]; i<(K_from[row-1]+K_sizes[row-1]); i++)
  {

    if (K_cols[i] == (col-1))
    {
      K_val[i] += val ; 
      return ;
    }

    if (K_cols[i] < 0)
    {
      K_cols[i] = (col-1) ;
      K_val[i] = val ; 
      return ;
    }
  }
  fprintf(stderr,"Addition of [%i,%i] to K failed (%e)\n",row,col,val);
  return; /* we should NOT reach this point */
}

/** Decomposition to L/U
 * @param a matrix (will be modified!)
 * @param index index vector
 * @param d modified index status (-1/+1)
 * @return status
 */
int femLUdecomp(a, index, size)
  double **a ;
  int     *index ;
  int      size ;
{
	int rv = 0;
	long i, j, k;
  long imax = 0  ;
	long n;
	double big,dum,sum,temp;
	double *vv ;

  n = size ;
	vv = NULL ;
  if ((vv= (double *)malloc(size*sizeof(double))) == NULL) goto memFree;
  for (i=0; i<size; i++) vv[i] = 0.0 ;

	for (i=0; i<n; i++)
	{
		big = 0.0 ;

		for (j=0; j<n; j++)
		{
			if ((temp=fabs((a[i][j]))) > big) 
      {
        big = temp;
      }
		}
		if (big == 0.0)
		{
			/* singular matrix */
			rv = -1 ; goto memFree;
		}
		vv[i] = (1.0/big) ;
	}

	for (j=0; j<n; j++)
	{
		for (i=0; i<j; i++)
		{
			sum = a[i][j];
			for (k=0; k<i; k++)
			{
				sum -= ( a[i][k]*a[k][j] ) ;
			}
			a[i][j] = sum ;
		}

		big = 0.0 ;
		for (i=j; i<n; i++)
		{
			sum = a[i][j];
			for (k=0; k<j; k++)
			{
				sum -= ( a[i][k]*a[k][j] ) ;
			}
			a[i][j] =  sum ;

			if ((dum=(vv[i]*fabs(sum))) >= big)
			{
				big  = dum ;
				imax = i ;
			}
		}

		if (j != imax) /* !? */
		{
			for (k=0; k<n; k++)
			{
				dum = a[imax][k] ;
				a[imax][k] = a[j][k] ;
				a[j][k] = dum ;
			}
			vv[imax] = vv[j] ;
		}

		index[j] =  imax ;

		if (a[j][j] == 0.0) {a[j][j] = 1e-10 ;}

		if (j != (n-1)) /* !? */
		{
			dum = 1.0 / a[j][j] ;
			for (i = (j+1); i<n; i++)
			{
				a[i][j] *= dum ;
			}
		}
	}

memFree:
	free(vv); vv = NULL ;
	return(rv);
}


/** Decomposition to L/U
 * @param a matrix (will be modified!)
 * @param index index vector
 * @param b right hand side/result vector (will be modified!)
 * @return status
 */
int femLUback(a, index, b, size)
  double **a ;
  int     *index ;
  double  *b ;
  int      size ;
{
	int    rv = 0 ;
	long   i,ii,ip,j;
	long   n ;
	double sum ;

	ii = -1 ;
  n = size ;

	for (i=0; i<n; i++)
	{
		ip  = index[i];
		sum = b[ip];
		b[ip] = b[i];

		if (ii > -1) /* means ii > 0 */
		{
			for (j=ii; j<=i-1; j++)
			{
				sum -= a[i][j]*b[j] ;
			}
		}
		else
		{
			if (sum)
			{
				ii = i ;
			}
		}
		b[i] = sum;
	}

	for (i=(n-1); i>=0; i--) 
	{
		sum = b[i];
		for (j=i+1; j<n; j++)
		{
			sum -= a[i][j]*b[j] ;
		}
		b[i] = (sum / a[i][i] ) ;
	}

	return(rv);
}


/** Inversion of "a" matrix using L/U
 * @param a matrix (will be modified!)
 * @return status
 */
int femLUinverse(a, size)
  double **a ;
  int      size;
{
	int       rv = 0 ;
  long      i,j ;
	long      n ;
	double  *col ;
	int     *index ;
	double  **b ;

	col = NULL  ;
	index = NULL;
	b = NULL ;
  n = size ;

  if ((col= (double *)malloc(size*sizeof(double))) == NULL) goto memFree;
  if ((index= (int *)malloc(size*sizeof(int))) == NULL) goto memFree;
  if ((b = (double **)malloc(size*sizeof(double *))) != NULL)
  {
    for (i=0; i<size; i++) b[i] = (double *)malloc(size*sizeof(double)) ;
  } else { exit(-1); }

  for (i=0; i<size; i++)
  {
    col[i] = 0.0 ;
    index[i] = 0.0 ;
    for (j=0; j<size; j++) b[i][j] = 0.0 ;
  }


	if ((rv = femLUdecomp(a, index, size)) != 0){goto memFree;}

  for (j=0; j<n; j++)
  {
    for (i=0; i<n; i++) { col[i] = 0.0 ; }

    col[j] = 1.0 ;
	  if ((rv = femLUback(a, index, col,size)) != 0){goto memFree;}

    for (i=0; i<n; i++) { b[i][j] = col[i] ; }
  }

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      a[i][j] = b[i][j]  ;
    }
  }

memFree:
	free(col); col = NULL ;
	free(index); index = NULL ;
    for (i=0; i<12; i++){free(b[i]); b[i]=NULL;} free(b); b=NULL;
	return(rv);
}


/** Support function for Ke computation */
double min_4(xi)
double *xi;
{
  int i; 
  double min ;
  
  min = xi[0] ;
  for (i=1; i<=4; i++) if (xi[i] < min) min = xi[i] ;
  return(min);
}

double max_4(xi)
double *xi;
{
  int i; 
  double max ;
  
  max = xi[0] ;
  for (i=1; i<=4; i++) if (xi[i] > max) max = xi[i] ;
  return(max);
}

/** Compute position of Fe member [1...3*n], node is from 1 (!) */
int pos_Fe(node, dir)
int node;
int dir;
{
  return((node-1)*3 + dir);
}


/* computes stiffness matrix of the structure */
void stiff(eg, lc)
int eg;
int lc;
{
  int i, j, k, ii, jj, m ;
  double xi[4];
  double yi[4];
  double A, tuh, xjj, yjj ;
  double a, b, c, d, tx, ty ;
  double xj[] = {-0.5774,  0.5774, 0.5774, -0.5774} ;
  double yj[] = {-0.5774, -0.5774, 0.5774,  0.5774} ;
  int   nodes[4] ;
  int   loc[12];

  for (i=0; i<n_elems; i++)
  {
    if ((e_g[i] > 0)&&(e_g[i] != eg)) continue;
    for (m=0; m<12; m++) 
    {
      fe[m] = 0.0 ;
      for (k=0;k<12;k++)
      {
        S[m][k] = 0.0 ;
        ke[m][k] = 0.0 ;
        SBDB[m][k] = 0.0 ;
        SBDBS[m][k] = 0.0 ;
      }
      for (k=0;k<3;k++)
      {
        SB[m][k] = 0.0 ;
        SBD[m][k] = 0.0 ;
      }
    }
    nodes[0] = n1[i]; nodes[1] = n2[i]; nodes[2] = n3[i]; nodes[3] = n4[i]; 
    for (ii=0; ii<4; ii++)
    {
      xi[ii] = x_i[nodes[ii]-1] ;
      yi[ii] = y_i[nodes[ii]-1] ;
    }

    A=(xi[1]-xi[0])*(yi[2]-yi[1]);
    if (A <= 0.0) {continue;} /* oops, zero area element */

	  a = min_4(xi);
	  b = max_4(xi);
	  c = min_4(yi);
	  d = max_4(yi);
    
    tuh = (double)((E[i]*pow(tl[i],3))/(12.0*(1.0-pow(nu[i],2)))) ;

    D[0][0] = tuh ;
    D[0][1] = tuh * nu[i] ;
    D[1][0] = D[0][1] ;
    D[1][1] = tuh ;
    D[2][2] = tuh*0.5*(1-nu[i]) ;

    for (j=0; j<4; j++)
    {
      S[j*3][0]=1.0;
      S[j*3][1]=xi[j];
      S[j*3][2]=yi[j];
      S[j*3][3]=pow(xi[j],2);
      S[j*3][4]=xi[j]*yi[j];
      S[j*3][5]=pow(yi[j],2);
      S[j*3][6]=pow(xi[j],3);
      S[j*3][7]=pow(xi[j],2)*yi[j];
      S[j*3][8]=xi[j]*pow(yi[j],2);
      S[j*3][9]=pow(yi[j],3);
      S[j*3][10]=pow(xi[j],3)*yi[j];
      S[j*3][11]=xi[j]*pow(yi[j],3);

      S[j*3+1][0]=0;
      S[j*3+1][1]=0;
      S[j*3+1][2]=1;
      S[j*3+1][3]=0;
      S[j*3+1][4]=xi[j];
      S[j*3+1][5]=2*yi[j];
      S[j*3+1][6]=0;
      S[j*3+1][7]=pow(xi[j],2);
      S[j*3+1][8]=2*xi[j]*yi[j];
      S[j*3+1][9]=3*pow(yi[j],2);
      S[j*3+1][10]=pow(xi[j],3);
      S[j*3+1][11]=3*xi[j]*pow(yi[j],2);

      S[j*3+2][0]=0;
      S[j*3+2][1]=-1;
      S[j*3+2][2]=0;
      S[j*3+2][3]=-2*xi[j];
      S[j*3+2][4]=-yi[j];
      S[j*3+2][5]=0;
      S[j*3+2][6]=-3*pow(xi[j],2);
      S[j*3+2][7]=-2*xi[j]*yi[j];
      S[j*3+2][8]=-pow(yi[j],2);
      S[j*3+2][9]=0;
      S[j*3+2][10]=-3*pow(xi[j],2)*yi[j];
      S[j*3+2][11]=-pow(yi[j],3);
    }

    femLUinverse(S,12) ; /* S inversion */

    /* integration start: */
    for (jj=0; jj<4; jj++)
    {
	    tx = (2.0*xj[jj]-a-b)/(b-a);
	    ty = (2.0*yj[jj]-c-d)/(d-c);

	    xjj = (((b-a)*tx)+b+a)/2.0;
	    yjj = (((d-c)*ty)+d+c)/2.0;

      B[0][3] = -2.0 ;
      B[0][6] = -6.0*xjj ;
      B[0][7] = -2.0*yjj ;
      B[0][10] = -6.0*xjj*yjj ;
      B[1][5] = -2.0 ;
      B[1][8] = -2.0*xjj ;
      B[1][9] = -6.0*yjj ;
      B[1][11] = -6.0*xjj*yjj ;
      B[2][4] = -2.0 ;
      B[2][7] = -4.0*xjj ;
      B[2][8] = -4.0*yjj ;
      B[2][10] = -6.0*xjj*xjj ;
      B[2][11] = -6.0*yjj*yjj ;

      /* ST * BT */
      tuh = 0.0 ;
      for (k=0; k<12; k++)
      {
        for (m=0; m<3; m++)
        {
          for (j=0; j<12; j++)
          {
            tuh += (S[j][k]*B[m][j]) ; /* transposition ! */
          }
          SB[k][m] = tuh ; tuh = 0.0 ;
        }
      }

      /* ST * BT * D */
      tuh = 0.0 ;
      for (k=0; k<12; k++) /* first dim */
      {
        for (m=0; m<3; m++) /* last dim */
        {
          for (j=0; j<3; j++) /* middle dim */
          {
            tuh += SB[k][j]*D[j][m] ;
          }
          SBD[k][m] = tuh ; tuh = 0.0 ;
        }
      }

      /* ST * BT * D * B */
      tuh = 0.0 ;
      for (k=0; k<12; k++) /* first dim */
      {
        for (m=0; m<12; m++) /* last dim */
        {
          for (j=0; j<3; j++) /* middle dim */
          {
            tuh += SBD[k][j]*B[j][m] ;
          }
          SBDB[k][m] = tuh ; tuh = 0.0 ;
        }
      }

      /* ST * BT * D * B * S */
      tuh = 0.0 ;
      for (k=0; k<12; k++) /* first dim */
      {
        for (m=0; m<12; m++) /* last dim */
        {
          for (j=0; j<12; j++) /* middle dim */
          {
            tuh += SBDB[k][j]*S[j][m] ;
          }
          SBDBS[k][m] = tuh ; tuh = 0.0 ;
        }
      }

      /* addition to the Ke (ke): */
      for (k=0; k<12; k++)
      {
        for (m=0; m<12; m++)
        {
          ke[k][m] += SBDBS[k][m] ;
        }
      }
    } /* end of integration */
    /* Multiplication bu area */
    for (k=0; k<12; k++)
    {
      for (m=0; m<12; m++)
      {
        ke[k][m] *= A ;
      }
    }

    /* TODO: print K matrix !!! */
    printf("\nKe[%i] in local coordinates:\n",i+1);
    for (k=0; k<12; k++) { for (m=0; m<12; m++) { printf("%10.3e ",ke[k][m]); } printf("\n"); }


    for (k=0; k<4; k++) /* localisation field */
    {
      for (m=0; m<3; m++)
      {
        loc[k*3+m] = pos_Fe(nodes[k], m+1) ;
      }
    }
    for (k=0; k<12; k++) /* actual localisation */
    {
      for (m=0; m<12; m++)
      {
        md_K_add(loc[k], loc[m], ke[k][m]) ;
      }
    }
  }
}

/* supports */
void add_one_disp(node, dir, val)
int   node;
int   dir;
double val;
{
  int i,j,n ;
  int row ;

  row = (node-1)*3 + dir - 1 ;  /* note: -1 ? */
  n = K_size ;

  for (i=0; i<n; i++)
  {
    for (j=K_from[i]; j<(K_from[i]+K_sizes[i]); j++)
    {
      if (K_cols[i] == (row))
      {
        if (K_cols[i] < 0) {break;}
        F_val[K_cols[i]] += K_val[i]*val ;
        K_val[i] = 0.0 ;
        break ;
      }
    }
  }

  for (i=K_from[row]; i<K_from[row]+K_sizes[row]; i++)
  {
    if (K_cols[i] < 0) {break;}
    if (K_cols[i] == row) 
    {
      K_val[i] = 1.0 ;
    }
    else
    {
      K_val[i] = 0.0 ;
    }
  }
  u_val[row] = val ;
  F_val[row] = val ; /* it will destroy any force in this place */
}


/* forces */
void add_one_force(node, dir, val)
int   node;
int   dir;
double val;
{
  int row ;

  row = (node-1)*3 + dir - 1 ;  /* note: -1 ? */

  F_val[row] += val ; 
}

void disps_and_loads(lc)
int lc;
{
  int i ;

  for (i=0; i<n_nfors; i++)
  {
    if ((f_g[i] != lc)) continue; /* only lc-th step is considered */

    add_one_force(f_n[i], f_d[i], f_v[i]);
  }

  for (i=0; i<n_disps; i++)
  {
    if ((d_g[i] > 0)&&(d_g[i] != lc)) continue;
    add_one_disp(d_n[i], d_d[i], d_v[i]);
  }
}

void print_f_f(void)
{
  int j ;

  printf("\nK |F:\n");
  for (j=0; j<(3*n_nodes); j++) { printf("| %10.3f\n", F_val[j]); }
}

/** Frees all allocated data */
void free_data()
{
  int i ;

  if (S != NULL)
  {
    for (i=0; i<12; i++){free(S[i]); S[i]=NULL;} free(S); S=NULL;
  }
  if (n_nodes > 0)
  {
    free(x_i); free(y_i);
  }
  if (n_elems > 0)
  {
    free(x_i); free(y_i); free(n1); free(n2); free(n3); free(n4);free(E);free(nu);free(tl);free(e_g);
  }
  if (n_disps > 0)
  {
    free(d_n);free(d_d);free(d_v);free(d_g);
  }
  if (n_nfors > 0)
  {
    free(f_n);free(f_d);free(f_v);free(f_g);
  }
}


void u_to_ue(s,c)
double s ;
double c ;
{
  int i,j ;
  double fval ;

  for (i=0; i<6; i++) { ue[i] = 0.0 ; }
  
  for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      fval += ue[j] ;
    }
    ue[i] = fval ;
  }
}


/** computes and prints structure description */
void res_struct(fw)
FILE *fw;
{
  int i;

  fprintf(fw,"\nDSLAB FINAL REPORT\n");

  fprintf(fw,"\nNodes:\n");
  fprintf(fw," Num     X            Y:\n");
  for (i=0; i<n_nodes; i++)
  {
    fprintf(fw," %3d %e %e\n",i+1, x_i[i], y_i[i]);
  }

  fprintf(fw,"\nElements:\n");
  fprintf(fw," Num Node1 Node2 Node3 Node4      E             A          nu:           t:           Grp:\n");
  for (i=0; i<n_elems; i++)
  {
    fprintf(fw," %3d   %3d   %3d   %3d   %3d %e %e %e ",i+1,n1[i],n2[i],n3[i],n4[i],E[i],nu[i],tl[i]);
    fprintf(fw," %2d\n",e_g[i]);
  }

  fprintf(fw,"\nSupports in nodes:\n");
  fprintf(fw," Num Node Dir Size:        LC:\n");
  for (i=0; i<n_disps; i++)
  {
    fprintf(fw," %3d  %3d ",i+1, d_n[i]);
    switch(d_d[i])
    {
      case 1: fprintf(fw," W ");break;
      case 2: fprintf(fw,"Fix");break;
      case 3: fprintf(fw,"Fiy");break;
      default: fprintf(fw,"unknown!");break;
    }
    fprintf(fw," %e %3d\n",d_v[i],d_g[i]);
  }

  fprintf(fw,"\nLoads in nodes:\n");
  fprintf(fw," Num Node Dir Size:         LC:\n");
  for (i=0; i<n_nfors; i++)
  {
    fprintf(fw," %3d  %3d ",i+1, f_n[i]);
    switch(f_d[i])
    {
      case 1: fprintf(fw,"\\|/");break;
      case 2: fprintf(fw,"Mx ");break;
      case 3: fprintf(fw,"My ");break;
      default: fprintf(fw,"unknown!");break;
    }
    fprintf(fw," %e %3d\n",f_v[i],f_g[i]);
  }

}

/** computes and prints results */
void results(fw)
FILE *fw;
{
  int i;

  fprintf(fw,"\nDEFORMATIONS:\n");
  fprintf(fw," Node  X            Y            Rotation:\n");
  for (i=0; i<n_nodes; i++)
  {
    fprintf(fw," %4d %e %e %e\n", i+1, u_val[3*i], u_val[3*i+1], u_val[3*i+2]);
  }

  for (i=0; i<n_disps; i++)
  {
    d_v[i] = 0.0 ; /* going to reuse it for reactions! */
  }

#if 0
  fprintf(fw,"\nINTERNAL FORCES (beginning and end of elements):\n");
  fprintf(fw," Element  FX1      FY1       M1        FX2        FY2       M2:\n");
  for (i=0; i<n_elems; i++)
  {
    /*res_loc(i);*/
    fprintf(fw," %3d %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e\n",
         i+1, fe[0], fe[1], fe[2], fe[3], fe[4], fe[5] );
  }

  fprintf(fw,"\nREACTIONS:\n");

  fprintf(fw," Num Node Dir Size:\n");
  for (i=0; i<n_disps; i++)
  {
    fprintf(fw," %3d  %3d ",i+1, d_n[i]);
    switch(d_d[i])
    {
      case 1: fprintf(fw,"-->");break;
      case 2: fprintf(fw,"/|\\");break;
      case 3: fprintf(fw,"__^");break;
      default: fprintf(fw,"unknown!");break;
    }
    fprintf(fw," %e\n",(-1.0)*d_v[i]);
  }
#endif
}



#ifndef GR
int main(argc, argv)
int argc ;
char *argv[];
{
  FILE *fw = NULL ;
  FILE *fo = NULL ;
#ifdef UI
	char  name[16] ;
	int   i ;
#endif

  fprintf(stderr,"\nDSLAB 0.0.1: slab solver for statics.\n");
  fprintf(stderr,"  See for details: http://github.com/jurabr/ddfor\n\n");

  if (argc < 2)
  {
#ifdef UI
    /* standard input */
    fprintf(stderr,"\nInteractive input:\n\n");
    if (read_data(stdin) != 0) 
    {
      fprintf(stderr,"\nProgram terminated!\n");
      exit(-1);
    }
		else
		{
			for (i=0; i<16; i++){name[i] = '\0';}
			fprintf(stderr,"\nName of file to save input (0 to skip):\n");
			if (fscanf(stdin,"%15s", name) > 0)
			{
				if (strlen(name) > 0)
				{
					if ((name[0]!=' ')&&(name[0]!='\n')&&(name[0] !='0'))
					{
    				if ((fw=fopen(name,"w")) != NULL)
						{
							write_idata(fw);
							fclose(fw);
						}
					}
				}
			}
		}
#else
    fprintf(stderr,"No input data file!\n");
    return(-1);
#endif
  }
  else
  {
    /* read from file */
    if ((fw=fopen(argv[1],"r")) == NULL)
    {
      fprintf(stderr,"Failed to open input file!\n");
      return(-1);
    }
    else
    {
      if (read_data(fw) != 0)
      {
        fprintf(stderr,"\nData input failed. Program terminated!\n");
        exit(-1);
      }
      else
      {
        fclose(fw);
      }
    }
  }

  if (argc < 3)
  {
    /* standard output */
    fo = stdout;
  }
  else
  {
    /* write to file */
    if ((argv[2][0] != '0')&&(argv[2][0] != '-'))
    {
      if ((fo=fopen(argv[2],"w")) == NULL)
      {
        fprintf(stderr,"Failed to open output file!\n");
        fo = stdout ;
      }
    } else {
      fo = NULL ;
    }
  }


  if (alloc_kf() != 0 )
  {
    free_data();
    return(-1);
  }

  	stiff(0,0); 
  	disps_and_loads(0);

  	fprintf(stderr,"\nSolution: \n");
  	solve_eqs();
  	fprintf(stderr,"End of solution. \n");

  	if (fo != NULL) 
  	{ 
    	res_struct(fo);
    	results(fo);
    	fclose(fo);
  	}
  
  free_sol_data();
  free_data();
  return(0);
}
#endif

/* end of ddfor.c */
