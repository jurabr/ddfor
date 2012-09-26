/*
 * File name: ddfor.c
 * Date:      2012/07/17 20:10
 * Author:    Jiri Brozovsky

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
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
#endif
#include <math.h>

int    n_nodes = 0 ;
int    n_elems = 0 ;
int    n_disps = 0 ;
int    n_nfors = 0 ;
int    n_eload = 0 ;

float *x_i = NULL ;
float *y_i = NULL ;

float *E = NULL ;
float *A = NULL ;
float *I = NULL ;
int   *n1 = NULL ;
int   *n2 = NULL ;
int   *type = NULL ;

int   *d_n = NULL ; /* node */
int   *d_d = NULL ; /* direction 1=x 2=y 3=rot */
float *d_v = NULL ; /* size */

int   *f_n = NULL ; /* node */
int   *f_d = NULL ; /* direction 1=fx 2=fy 3=m */
float *f_v = NULL ; /* size */

int   *l_e = NULL ; /* node */
int   *l_d = NULL ; /* direction 1=x 2=y, 3=x global, 4=y global */
float *l_v1 = NULL ; /* size at beginning */
float *l_v2 = NULL ; /* size at end */

/* solution variables: */
double  ke[6][6] ;
double  keg[6][6] ;
double  T[6][6] ;
double  fe[6];
double  feg[6];
double  ueg[6];
double  ue[6];

int    K_len    = 0 ;
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
  fprintf(stderr,"  Declared %d nodes.\n", n_nodes);

  /* allocate data for nodes */
  if ((x_i=(float *)malloc(n_nodes*sizeof(float))) == NULL)
  {
    fprintf(stderr,"Can't allocate X coordinates!\n");
    return(-2);
  }

  fprintf(stderr,"Coordinates of nodes (x y):\n");
  if ((y_i=(float *)malloc(n_nodes*sizeof(float))) == NULL)
  {
    free(x_i);
    fprintf(stderr,"Can't allocate Y coordinates!\n");
    return(-2);
  }

  /* read coordinates of nodes */
  for (i=0; i<n_nodes; i++)
  {
    fscanf(fw,"%e %e", &x_i[i], &y_i[i]) ;
    fprintf(stderr," %d %e %e\n", i+1, x_i[i], y_i[i]) ;
  }
  fprintf(stderr,"  Have %d coordinates.\n",n_nodes);


  /* number of elements */
  fprintf(stderr,"Number of elements:\n");
  fscanf(fw,"%d", &n_elems);
  if (n_elems <= 0)
  {
    free(x_i); free(y_i); 
    fprintf(stderr,"No elements!\n");
    return(-1);
  }
  fprintf(stderr,"  Declared %d elements.\n", n_elems);

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
  if ((E=(float *)malloc(n_elems*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);
    fprintf(stderr,"Can't allocate stifnesses!\n");
    return(-2);
  }
  if ((A=(float *)malloc(n_elems*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);
    fprintf(stderr,"Can't allocate areas!\n");
    return(-2);
  }
  if ((I=(float *)malloc(n_elems*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);
    fprintf(stderr,"Can't allocate moments of inertia!\n");
    return(-2);
  }
  if ((type=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    fprintf(stderr,"Can't allocate types!\n");
    return(-2);
  }

  fprintf(stderr,"Element data (type node1 node2 E A I):\n");
  /* read element data */
  for (i=0; i<n_elems; i++)
  {
    fscanf(fw,"%d %d %d %e %e %e",&type[i], &n1[i], &n2[i], &E[i], &A[i], &I[i]) ;
    fprintf(stderr,"%d %d %d %e %e %e\n",type[i], n1[i], n2[i], E[i], A[i], I[i]) ;
  }
  fprintf(stderr,"  Have %d elements.\n",n_elems);

  /* supports */
  fprintf(stderr,"Number of supports:\n");
  fscanf(fw,"%d", &n_disps);
  if (n_disps <= 0)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    fprintf(stderr,"No supports!\n");
    return(-1);
  }
  fprintf(stderr,"  Declared %d supports.\n", n_disps);

  if ((d_n=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    fprintf(stderr,"Can't allocate nodes for supports!\n");
    return(-2);
  }
  if ((d_d=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);
    fprintf(stderr,"Can't allocate types of supports!\n");
    return(-2);
  }
  if ((d_v=(float *)malloc(n_disps*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);
    fprintf(stderr,"Can't allocate sizes of supports!\n");
    return(-2);
  }
  
  fprintf(stderr,"Supports data (node direction size):\n");
  /* read supports data */
  for (i=0; i<n_disps; i++)
  {
    fscanf(fw,"%d %d %e",&d_n[i], &d_d[i], &d_v[i]) ;
    fprintf(stderr,"%d %d %e\n",d_n[i], d_d[i], d_v[i]) ;
  }
  fprintf(stderr,"  Have %d supports.\n",n_disps);


  /* forces in nodes */
  fprintf(stderr,"Number of forces in nodes:\n");
  fscanf(fw,"%d", &n_nfors);
  if (n_nfors <= 0)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    fprintf(stderr,"No forces in nodes!\n");
    return(-1);
  }
  fprintf(stderr,"  Declared %d forces in nodes.\n", n_nfors);

  if ((f_n=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    fprintf(stderr,"Can't allocate nodes for forces!\n");
    return(-2);
  }
  if ((f_d=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    free(f_n);
    fprintf(stderr,"Can't allocate types of forces!\n");
    return(-2);
  }
  if ((f_v=(float *)malloc(n_nfors*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    free(f_n);free(f_d);
    fprintf(stderr,"Can't allocate sizes of forces!\n");
    return(-2);
  }
  
  fprintf(stderr,"Forces in nodes data (node direction size):\n");
  /* read supports data */
  for (i=0; i<n_nfors; i++)
  {
    fscanf(fw,"%d %d %e",&f_n[i], &f_d[i], &f_v[i]) ;
    fprintf(stderr,"%d %d %e\n",f_n[i], f_d[i], f_v[i]) ;
  }
  fprintf(stderr,"  Have %d forces in nodes.\n",n_nfors);


  /* loads on elements */
  fprintf(stderr,"Number of element loads:\n");
  fscanf(fw,"%d", &n_eload);
  if (n_eload <= 0)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    free(f_n);free(f_d);free(f_v);
    fprintf(stderr,"No element loads!\n");
    return(-1);
  }
  fprintf(stderr,"  Declared %d element loads.\n", n_eload);

  if ((l_e=(int *)malloc(n_eload*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    free(f_n);free(f_d);free(f_v);
    fprintf(stderr,"Can't allocate elements for loads!\n");
    return(-2);
  }
  if ((l_d=(int *)malloc(n_eload*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);}
    free(l_e);
    fprintf(stderr,"Can't allocate types of loads!\n");
    return(-2);
  }
  if ((l_v1=(float *)malloc(n_eload*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);}
    free(l_e);free(l_d);
    fprintf(stderr,"Can't allocate sizes of loads!\n");
    return(-2);
  }
  if ((l_v2=(float *)malloc(n_eload*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);

    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);}
    free(l_e);free(l_d);free(l_v1);
    fprintf(stderr,"Can't allocate sizes of loads!\n");
    return(-2);
  }
  
  fprintf(stderr,"Element loads (element direction startsize endsize):\n");
  /* read supports data */
  for (i=0; i<n_eload; i++)
  {
    fscanf(fw,"%d %d %e %e",&l_e[i], &l_d[i], &l_v1[i], &l_v2[i]) ;
    fprintf(stderr,"%d %d %e %e\n",l_e[i], l_d[i], l_v1[i], l_v2[i]) ;
  }
  fprintf(stderr,"  Have %d element loads.\n",n_eload);

  return(0);
}

/** free K,F */
void free_sol_data()
{
  if (F_val != NULL )free(F_val);
  if (u_val != NULL )free(u_val);
  if (K_sizes != NULL )free(K_sizes);
  if (K_from != NULL )free(K_from);
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

  if ((K_sizes = (int *)malloc(3*n_nodes*sizeof(int)))   == NULL) { goto memFree ; }
  if ((K_from  = (int *)malloc(3*n_nodes*sizeof(int)))   == NULL) { goto memFree ; } 
  if ((F_val   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 
  if ((u_val   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 

  if ((M   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 
  if ((r   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 
  if ((z   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 
  if ((p   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 
  if ((q   = (double *)malloc(3*n_nodes*sizeof(double))) == NULL) { goto memFree ; } 

  K_len = 0 ;

  for (i=0; i<n_nodes*3; i++) 
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
				K_sizes[3*i+0]+=6; 
				K_sizes[3*i+1]+=6; 
				K_sizes[3*i+2]+=6; 
			}
      if ((n2[j]-1) == i) 
			{
				K_sizes[3*i+0]+=6; 
				K_sizes[3*i+1]+=6; 
				K_sizes[3*i+2]+=6; 
			}
    }
  }

  sum = 0 ;

  for (i=0; i<(n_nodes*3); i++)
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

  n = 3*n_nodes ;

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

/** Local stiffness matrix */
void stiff_loc(type, E, A, I, l)
int type ;
float E ;
float A ;
float I ;
double l ;
{
  int i, j ;
  float tuh ;

  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      ke[i][j] = 0.0 ;
      keg[i][j] = 0.0 ;
    }
  }

  /* common for all: */
  tuh = (double) ((E*A)/l) ;
  ke[0][0] = tuh ;
  ke[3][3] = tuh ;
  ke[0][3] =( -tuh) ;
  ke[3][0] =( -tuh) ;


  switch (type)
  {
    case 0: /* |--| */
      ke[1][1] = (12.0*E*I)/(l*l*l) ;
      ke[1][4] = (-12.0*E*I)/(l*l*l) ;
      ke[4][1] = (-12.0*E*I)/(l*l*l) ;
      ke[4][4] = (12.0*E*I)/(l*l*l) ;

      ke[1][2] = (-6.0*E*I)/(l*l) ;
      ke[1][5] = (-6.0*E*I)/(l*l) ;
      ke[2][1] = (-6.0*E*I)/(l*l) ;
      ke[2][4] = (6.0*E*I)/(l*l) ;

      ke[4][2] = (6.0*E*I)/(l*l) ;
      ke[4][5] = (6.0*E*I)/(l*l) ;
      ke[5][1] = (-6.0*E*I)/(l*l) ;
      ke[5][4] = (6.0*E*I)/(l*l) ;

      ke[2][2] = (4.0*E*I)/(l) ;
      ke[5][5] = (4.0*E*I)/(l) ;

      ke[2][5] = (2.0*E*I)/(l) ;
      ke[5][2] = (2.0*E*I)/(l) ;
      break;
    case 1: /* o--| */
      ke[1][1] = (3.0*E*I)/(l*l*l) ;
      ke[1][4] = (-3.0*E*I)/(l*l*l) ;
      ke[1][5] = (-3.0*E*I)/(l*l) ;

      ke[4][1] = (-3.0*E*I)/(l*l*l) ;
      ke[4][4] = (3.0*E*I)/(l*l*l) ;
      ke[4][5] = (3.0*E*I)/(l*l) ;

      ke[5][1] = (-3.0*E*I)/(l*l) ;
      ke[5][4] = (3.0*E*I)/(l*l) ;
      ke[5][5] = (3.0*E*I)/(l) ;
      break;
    case 2: /* |--o */ 
      ke[1][1] = (3.0*E*I)/(l*l*l) ;
      ke[1][2] = (-3.0*E*I)/(l*l) ;
      ke[1][4] = (-3.0*E*I)/(l*l*l) ;

      ke[2][1] = (-3.0*E*I)/(l*l) ;
      ke[2][2] = (3.0*E*I)/(l) ;
      ke[2][4] = (3.0*E*I)/(l*l) ;

      ke[4][1] = (-3.0*E*I)/(l*l*l) ;
      ke[4][2] = (3.0*E*I)/(l*l) ;
      ke[4][4] = (3.0*E*I)/(l*l*l) ;
      break;
    case 3: /* o--o .. nothing to do */
      break;
  }
}


/* set transformation matrix to zero */
void tran_zero()
{
  int i, j;
  for (i=0; i<6; i++) { for (j=0; j<6; j++) { T[i][j] = 0.0 ; } }
  T[2][2] = 1.0 ;
  T[5][5] = 1.0 ;
}

/** Transformation matrix */
void tran(s, c)
double s ;
double c ;
{
  tran_zero();
  T[0][0] = c ;
  T[0][1] = s ;
  T[1][0] = (-s) ;
  T[1][1] = c ;

  T[3][3] = c ;
  T[3][4] = s ;
  T[4][3] = (-s) ;
  T[4][4] = c ;
}

/* fils "ke" with content of "keg" */
void ke_switch()
{
  int i, j;
  for (i=0; i<6; i++) { for (j=0; j<6; j++) { ke[i][j] = keg[i][j] ; } }
}

void ke_to_keg(s, c)
double s ;
double c ;
{
  int i,j,k;
  float fval, kval ;

  tran(s, c);

  for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      fval += T[j][i] * fe[j] ;

      kval = 0.0 ;
      for (k=0; k<6; k++)
      {
        kval += T[k][i] * ke[k][j] ;
      }
      keg[i][j] = kval ;
    }
    feg[i] = fval ;
  }

  ke_switch();

  for (i=0; i<6; i++)
  {
    fval = 0.0 ;
    for (j=0; j<6; j++)
    {
      kval = 0.0 ;
      for (k=0; k<6; k++)
      {
        kval += ke[i][k] * T[k][j] ;
      }
      keg[i][j] = kval ;
    }
  }
}

/* puts data to the right place in K */
void md_K_add(row, col, val)
int row;
int col;
float val;
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

/* loads on nodes */
void one_eload(epos, na, nb, va, vb, L)
int epos;
double na;
double nb;
double va;
double vb;
double L;
{
	switch (type[epos])
  {
    case 0: /* |--| */
            fe[0]+=(-(2.0*na+1.0*nb)*L)/6.0 ;
            fe[1]+=(-(7.0*va+3.0*vb)*L)/20.0 ;
            fe[2]+=((3.0*va+2.0*vb)*L*L)/60.0 ;
            fe[3]+=(-(1.0*na+2.0*nb)*L)/6.0 ;
            fe[4]+=(-(3.0*va+7.0*vb)*L)/20.0 ;
            fe[5]+=(-(2.0*va+3.0*vb)*L*L)/60.0 ;
            break ;
    case 1: /* o--| */
            fe[0]+=(-(2.0*na+1.0*nb)*L)/6.0 ;
            fe[1]+=(-(4.0*va+11.0*vb)*L)/40.0 ;
            fe[3]+=(-(1.0*na+2.0*nb)*L)/6.0 ;
            fe[4]+=(-(16.0*va+9.0*vb)*L)/40.0 ;
            fe[5]+=((8.0*va+7.0*vb)*L*L)/120.0 ;
            break ;
    case 2: /* |--o */
            fe[0]+=(-(2.0*na+1.0*nb)*L)/6.0 ;
            fe[1]+=(-(16.0*va+9.0*vb)*L)/40.0 ;
            fe[3]+=(-(1.0*na+2.0*nb)*L)/6.0 ;
            fe[2]+=((8.0*va+7.0*vb)*L*L)/120.0 ;
            fe[4]+=(-(4.0*va+11.0*vb)*L)/40.0 ;
            break ;
    case 3: /* o--o */
            fe[0]+=(-(2.0*na+1.0*nb)*L)/6.0 ;
            fe[1]+=(-(7.0*va+3.0*vb)*L)/20.0 ;
            fe[3]+=(-(1.0*na+2.0*nb)*L)/6.0 ;
            fe[4]+=(-(3.0*va+7.0*vb)*L)/20.0 ;
            break ;
    default: return; break;
  }
}

/* computes stiffness matrix of the structure */
void stiff()
{
  int i, j, k, ii, jj, m ;
  float x1,y1, x2,y2, l, s, c ;


  for (i=0; i<n_elems; i++)
  {
	  for (m=0; m<6; m++) {fe[m] = 0.0 ; feg[m]=0.0;}
    x1 = x_i[n1[i]-1] ;
    y1 = y_i[n1[i]-1] ;
    x2 = x_i[n2[i]-1] ;
    y2 = y_i[n2[i]-1] ;
    
    l = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) ;
    if (l <= 0.0) {continue;} /* oops, zero length element */
    s = (y2-y1)/l ;
    c = (x2-x1)/l ;
    
    tran_zero();
    stiff_loc(type[i], E[i], A[i], I[i], (double)l) ;

		/* loads on elements: */
		for (m=0; m<n_elems; m++)
		{
			if ((l_e[m]-1) == i)
			{
				switch (l_d[m])
				{
					case 1: 
						one_eload(i, (double)l_v1[m], (double)l_v2[m], 0.0, 0.0, (double)l);
						break;
					case 2: 
						one_eload(i, 0.0, 0.0, (double)l_v1[m], (double)l_v2[m], (double)l);
						break;
				}
			}
		}

		/* ke, fe transformation: */
    ke_to_keg(s, c) ;

    /* localisation */
    for (k=0; k<6; k++)
    {
      if (k <3) { ii = n1[i]*3 + k - 2 ; }
      else      { ii = n2[i]*3 + k - 5 ; }

      F_val[ii-1] += feg[k] ;
  
      for (j=0; j<6; j++)
      {
        if (j <3) { jj = n1[i]*3 + j - 2 ; }
        else      { jj = n2[i]*3 + j - 5 ; }
  
        md_K_add(ii, jj, keg[k][j]) ;
      }
    }
  }
}

/* supports */
void add_one_disp(node, dir, val)
int   node;
int   dir;
float val;
{
  int i,j,n ;
  int row ;

  row = (node-1)*3 + dir - 1 ;  /* note: -1 ? */
  n = 3*n_nodes ;

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
float val;
{
  int row ;

  row = (node-1)*3 + dir - 1 ;  /* note: -1 ? */

  F_val[row] += val ; 
}

void disps_and_loads()
{
  int i ;

  for (i=0; i<n_nfors; i++)
  {
		add_one_force(f_n[i], f_d[i], f_v[i]);
  }

  for (i=0; i<n_disps; i++)
  {
		add_one_disp(d_n[i], d_d[i], d_v[i]);
  }
}

/** Frees all allocated data */
void free_data()
{
  if (n_nodes > 0)
  {
    free(x_i); free(y_i);
  }
  if (n_elems > 0)
  {
    free(n1); free(n2);free(E);free(A);free(I);
  }
  if (n_disps > 0)
  {
    free(d_n);free(d_d);free(d_v);
  }
  if (n_nfors > 0)
  {
    free(f_n);free(f_d);free(f_v);
  }
  if (n_eload > 0)
  {
    free(l_e);free(l_d);free(l_v1);free(l_v2);
  }
}

void u_to_ue(s,c)
double s ;
double c ;
{
	int i,j ;
	double fval ;

  tran(s, c);
	for (i=0; i<6; i++) { ue[i] = 0.0 ; }
  
	for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      fval += T[i][j] * ueg[j] ;
    }
    ue[i] = fval ;
  }
}

/** Find all element load in give direction */
void get_eloads(epos, dir, na, nb)
int epos;
int dir;
double *na;
double *nb;
{
  int i;
  
  *na = 0.0 ;
  *nb = 0.0 ;

  for (i=0; i<n_eload; i++)
  {
    if ( (l_e[i] == (epos+1)) && (dir == l_d[i]) )
    {
      *na += l_v1[i] ;
      *nb += l_v2[i] ;
    }
  }
}


/** Compute internal force (N, V, M) for given point of beam */
double in_force(type, epos, div, ppos)
int type; /* force type: 1=N, 2=V, 3=M */
int epos; /* element position */
int div;  /* number of divisions */
int ppos; /* number of computed point (0...div)*/
{
  double x1,x2,y1,y2 ;
  double Na,Nb, Va,Vb, Ma,Mb, L, lenx, lenxx, na, nb, no, nt ;
  double Xo = 0.0 ;

  Na = fe[0];
  Va = fe[1];
  Ma = fe[2];
  Nb = fe[3];
  Vb = fe[4];
  Mb = fe[5];

	x1 = x_i[n1[epos]-1] ;
  y1 = y_i[n1[epos]-1] ;
  x2 = x_i[n2[epos]-1] ;
  y2 = y_i[n2[epos]-1] ;
 

  L = sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) ) ;
  lenx  = L*((double)((double)ppos/(double)(div))) ;
  lenxx = L - lenx ;

  switch (type)
  {
    case 0 : return(lenx); break;
    case 1 : get_eloads(epos, 1, &na, &nb);
             no = na ; nt = nb - no ;
             Xo = lenx*no + 0.5*lenx*((nt*lenx/L)) ;
             return (Xo - (Na) ) ;
             break ;
    case 2 : get_eloads(epos, 2, &na, &nb);
             no = na ; nt = nb - no ;

             Xo = (no*L)/2 - (no*lenx) 
							 + ((nt*L)/6 - ((nt*lenx*lenx)/(2*L)) ) ;
             return (Xo  - ((Mb + Ma)/L)  ) ;
             break ;
    case 3 : get_eloads(epos, 2, &na, &nb);
             no = na ; nt = nb - no ;

             Xo =  (no*L*lenx)/2 - (no*lenx*lenx)/2 
                  + ((nt*L*lenx)/6.0 - (nt*lenx*lenx*lenx)/(6*L) ) ;
             return ((1.0)*(Xo + ((-Ma*lenxx+Mb*lenx)/L)  ) ) ;
             break ;
  }

  return(0.0);
}

/** Local results */
void res_loc(epos)
int epos ;
{
  int i, j, m ;
  float x1,y1, x2,y2, l, s, c, fval ;

	for (i=0; i<6; i++) 
	{ 
		fe[i]  = 0.0 ; 
		feg[i] = 0.0 ; 
		ue[i]  = 0.0 ; 
	}
    
	/* get initial stuff */
	for (i=0; i<3; i++)
	{
		ueg[i]   = u_val[(3*(n1[epos]-1))+i];
		ueg[i+3] = u_val[(3*(n2[epos]-1))+i];
	}

	x1 = x_i[n1[epos]-1] ;
  y1 = y_i[n1[epos]-1] ;
  x2 = x_i[n2[epos]-1] ;
  y2 = y_i[n2[epos]-1] ;
    
  l = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) ;
  if (l <= 0.0) {return;} /* oops, zero length element */
  s = (y2-y1)/l ;
  c = (x2-x1)/l ;

  u_to_ue(s,c);
  stiff_loc(type[epos], E[epos], A[epos], I[epos], (double)l) ;

	/* compute forces in nodes: */
	for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      fval += ke[i][j] * ue[j] ;
    }
    fe[i] -= fval ;
  }

	/* primary forces: */
	for (m=0; m<n_elems; m++)
	{
		if ((l_e[m]-1) == epos)
		{
			switch (l_d[m])
			{
				case 1: 
					one_eload(epos, (double)l_v1[m], (double)l_v2[m], 0.0, 0.0, (double)l);
					break;
				case 2: 
					one_eload(epos, 0.0, 0.0, (double)l_v1[m], (double)l_v2[m], (double)l);
					break;
			}
		}
	}

	/* get reactions: */
	for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      fval += T[j][i] * fe[j] ;
    }
    feg[i] += fval ;
  }

	for (i=0; i<n_disps; i++)
	{
		if (d_n[i] == n1[epos]) { d_v[i] +=  feg[d_d[i]-1] ; }
		if (d_n[i] == n2[epos]) { d_v[i] +=  feg[d_d[i]+2] ; }

		for (j=0; j<n_nfors; j++)
		{
			if ((d_n[i] == f_n[j]) && (d_d[i] == f_d[j]))
			{
				d_v[i] += f_v[j] ;
			}
		}
	}
}


/** computes and prints results */
void results(fw)
FILE *fw;
{
	int i;

	fprintf(fw,"\nFINAL REPORT\n");

	fprintf(fw,"\nNodes:\n");
	fprintf(fw," Num     X            Y:\n");
	for (i=0; i<n_nodes; i++)
	{
		fprintf(fw," %3d %e %e\n",i+1, x_i[i], y_i[i]);
	}

	fprintf(fw,"\nElements:\n");
	fprintf(fw," Num Node1 Node2     E             A           I:\n");
	for (i=0; i<n_elems; i++)
	{
		fprintf(fw," %3d   %3d   %3d %e %e %e ",i+1,n1[i],n2[i],E[i],A[i],I[i]);
		switch(type[i])
		{
			case 0: fprintf(fw,"|--|\n");break;
			case 1: fprintf(fw,"o--|\n");break;
			case 2: fprintf(fw,"|--o\n");break;
			case 3: fprintf(fw,"o--o\n");break;
			default: fprintf(fw,"unknown!\n");break;
		}
	}

	fprintf(fw,"\nSupports in nodes:\n");
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
		fprintf(fw," %e\n",d_v[i]);
	}

	fprintf(fw,"\nLoads in nodes:\n");
	fprintf(fw," Num Node Dir Size:\n");
	for (i=0; i<n_nfors; i++)
	{
		fprintf(fw," %3d  %3d ",i+1, f_n[i]);
		switch(f_d[i])
		{
			case 1: fprintf(fw,"-->");break;
			case 2: fprintf(fw,"/|\\");break;
			case 3: fprintf(fw,"__^");break;
			default: fprintf(fw,"unknown!");break;
		}
		fprintf(fw," %e\n",f_v[i]);
	}

	fprintf(fw,"\nElement loads:\n");
	fprintf(fw," Num Node Dir Sizes (start..end):\n");
	for (i=0; i<n_eload; i++)
	{
		fprintf(fw," %3d  %3d ",i+1, l_e[i]);
		switch(l_d[i])
		{
			case 1: fprintf(fw,">>>");break;
			case 2: fprintf(fw,"vvv");break;
			default: fprintf(fw,"unknown!");break;
		}
		fprintf(fw," %e..%e\n",l_v1[i],l_v2[i]);
	}

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

	fprintf(fw,"\nINTERNAL FORCES (beginning and end of elements):\n");
	fprintf(fw," Element  FX1      FY1       M1        FX2        FY2       M2:\n");
	for (i=0; i<n_elems; i++)
	{
		res_loc(i);
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
}

/** computes and prints internal forces in elements */
void eint_results(fw)
FILE *fw;
{
  int i,j ;
  int div = 10 ;

  fprintf(fw,"\n# DETAILED INTERNAL FORCES: #\n");
	fprintf(fw,"# Element  x      FX        FY        M\n");
	for (i=0; i<n_elems; i++)
	{
		res_loc(i);
    for (j=0; j<=div; j++)
    {
		fprintf(fw," %3d %2.3e %2.3e %2.3e %2.3e\n",
			 	i+1,
        in_force(0, i, div, j),
        in_force(1, i, div, j),
        in_force(2, i, div, j),
        in_force(3, i, div, j)
        );
    }
	  
    fprintf(fw,"\n");
	}
}

/** Compute internal force (N, V, M) for given point of beam
 *  This routine is meant for plotting tools */
void in_gfx(type, epos, div, ppos, mult, vx, vy)
int type; /* force type: 1=N, 2=V, 3=M */
int epos; /* element position */
int div;  /* number of divisions */
int ppos; /* number of computed point (0...div)*/
double mult;
double *vx;
double *vy;
{
  double x1,x2,y1,y2,c,s ;
  double Na,Nb, Va,Vb, Ma,Mb, L, lenx, lenxx, na, nb, no, nt ;
  double Xo = 0.0 ;
  double val = 0.0 ;

  Na = fe[0];
  Va = fe[1];
  Ma = fe[2];
  Nb = fe[3];
  Vb = fe[4];
  Mb = fe[5];

	x1 = x_i[n1[epos]-1] ;
  y1 = y_i[n1[epos]-1] ;
  x2 = x_i[n2[epos]-1] ;
  y2 = y_i[n2[epos]-1] ;
 

  L = sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) ) ;
  lenx  = L*((double)((double)ppos/(double)(div))) ;
  lenxx = L - lenx ;

  s = (y2-y1)/L ;
  c = (x2-x1)/L ;

  *vx = x1 + lenx*c ;
  *vy = y1 + lenx*s ;
  val = 0 ;

  switch (type)
  {
    case 0 : /* already computed */ break;
    case 1 : get_eloads(epos, 1, &na, &nb);
             no = na ; nt = nb - no ;
             Xo = lenx*no + 0.5*lenx*((nt*lenx/L)) ;
             val = (Xo - (Na) ) ;
             break ;
    case 2 : get_eloads(epos, 2, &na, &nb);
             no = na ; nt = nb - no ;

             Xo = (no*L)/2 - (no*lenx) 
							 + ((nt*L)/6 - ((nt*lenx*lenx)/(2*L)) ) ;
             val = (Xo  - ((Mb + Ma)/L)  ) ;
             break ;
    case 3 : get_eloads(epos, 2, &na, &nb);
             no = na ; nt = nb - no ;

             Xo =  (no*L*lenx)/2 - (no*lenx*lenx)/2 
                  + ((nt*L*lenx)/6.0 - (nt*lenx*lenx*lenx)/(6*L) ) ;
             val = ((-1.0)*(Xo + ((-Ma*lenxx+Mb*lenx)/L)  ) ) ;
             break ;
  }

  *vx = *vx - val*mult*s ;
  *vy = *vy + val*mult*c ;
}

/** prints internal forces in elements gfx-friendly form */
void gfx_results(fw)
FILE *fw;
{
  int i,j ;
  int div = 10 ;
  double x,y, nx,ny, vx,vy, mx,my ;
  double size = 0.0;
  double mult = 0.0;
  double dmult = 0.0;

  /* compute multiplier(s) first */
  for (i=0; i<n_nodes; i++)
  {
    if (fabs(x_i[i]) > size) {size = fabs(x_i[i]);}
    if (fabs(y_i[i]) > size) {size = fabs(y_i[i]);}
    for (j=(i*3); j<((i*3)+3); j++) /* unused just now */
    {
      if (fabs(u_val[j]) > dmult) {dmult = fabs(u_val[j]);}
    }
  }
  if (dmult > 1e-8) { dmult = 0.3*size/dmult; }

	for (i=0; i<n_disps; i++)
  {
    if (fabs(d_v[i]) > mult) {mult = fabs(d_v[i]);}
  }
  if (mult > 1e-6){mult = 0.3*size/mult ;}

  fprintf(fw,"\n# DETAILED INTERNAL FORCES (plot data): #\n");
	fprintf(fw,"# Element  x       y      Nx    Ny      Vx    Vy      Mx      My\n");
	for (i=0; i<n_elems; i++)
	{
		res_loc(i);
    for (j=0; j<=div; j++)
    {
      in_gfx(0, i, div, j, mult, &x,  &y);
      in_gfx(1, i, div, j, mult, &nx, &ny);
      in_gfx(2, i, div, j, mult, &vx, &vy);
      in_gfx(3, i, div, j, mult, &mx, &my);

      /* first line - for nicer plot... */
      if (j == 0)
      {
        fprintf(fw,
        " %3d %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e\n",
			 	i+1, x,y, x,y,x,y,x,y 
        );
      }

      /* actual data output */
		  fprintf(fw,
        " %3d %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e\n",
			 	i+1, x,y, nx,ny,vx,vy,mx,my 
      );
      
      /* last line - for nicer plot... */
      if (j==div)
      {
        fprintf(fw,
        " %3d %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e %2.3e\n",
			 	i+1, x,y, x,y,x,y,x,y 
        );
      }
    }
	  
    fprintf(fw,"\n");
	}
}

/** Find (approximation of) extreme values of N, V, M on beam */
int beam_max(type, epos, div, nmax, npos, vmax, vpos, mmax, mpos)
int type; /* force type: 1=N, 2=V, 3=M */
int epos; /* element position */
int div;  /* number of divisions */
double *nmax; /* maximal N */
double *npos; /* position of max N absolute from 1st node */
double *vmax; /* maximal V */
double *vpos; /* position of max V absolute from 1st node */
double *mmax; /* maximal M */
double *mpos; /* position of max M absolute from 1st node */
{
  int i ;

  *nmax = 0.0 ;
  *vmax = 0.0 ;
  *mmax = 0.0 ;

  for (i=0; i<=div; i++)
  {
    /* TODO */
  }

  /* TODO code here */
  return(0);
}

int main(argc, argv)
int argc ;
char *argv[];
{
  FILE *fw = NULL ;
  FILE *fo = NULL ;
  FILE *fd = NULL ;
  FILE *fp = NULL ;

	fprintf(stderr,"\nDDFOR 1.0.1: direct stiffness method solver for statics of 2D frames.\n");
	fprintf(stderr,"  See for details: http://github.com/jurabr/ddfor\n\n");

  if (argc < 2)
  {
#ifdef LARGE
    /* standard input */
    fprintf(stderr,"\nInteractive input:\n\n");
    read_data(stdin);
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
      read_data(fw);
      fclose(fw);
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
    if ((fo=fopen(argv[2],"w")) == NULL)
    {
      fprintf(stderr,"Failed to open output file!\n");
			fo = stdout ;
    }
  }

  if (argc < 4)
  {
    /* no output */
    fd = NULL;
  }
  else
  {
    /* write to file */
    if ((fd=fopen(argv[3],"w")) == NULL)
    {
      fprintf(stderr,"Failed to open output file!\n");
			fd = NULL ;
    }
  }

  if (argc < 5)
  {
    /* no output */
    fp = NULL;
  }
  else
  {
    /* write to file */
    if ((fp=fopen(argv[4],"w")) == NULL)
    {
      fprintf(stderr,"Failed to open gfx file!\n");
			fp = NULL ;
    }
  }


  if (alloc_kf() != 0 )
  {
    free_data();
    return(-1);
  }

  stiff(); 
	disps_and_loads();

	fprintf(stderr,"\nSolution: \n");
  solve_eqs();
	fprintf(stderr,"End of solution. \n");

  results(fo);
	fclose(fo);

  if (fd != NULL) 
  { 
    eint_results(fd); 
    fclose(fd);
  }

  if (fp != NULL) 
  { 
    gfx_results(fp); 
    fclose(fp);
  }


  free_sol_data();
  free_data();
  return(0);
}

/* end of ddfor.c */
