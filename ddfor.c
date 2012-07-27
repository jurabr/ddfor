/*
 * File name: ddfor.c
 * Date:      2012/07/17 20:10
 * Author:    Jiri Brozovsky
 */

#include <stdio.h>
/*#include <stdlib.h>*/
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
float  ke[6][6] ;
float  keg[6][6] ;
float  T[6][6] ;
float  fe[6];
float  feg[6];
float  ueg[6];

int    K_len    = 0 ;
int   *K_sizes  = NULL ; /* lenghts of K's rows */
int   *K_from   = NULL ;
int   *K_cols   = NULL ;
float *K_val    = NULL ;
float *F_val    = NULL ;
float *u_val    = NULL ;

float *M    = NULL ;
float *r    = NULL ;
float *z    = NULL ;
float *p    = NULL ;
float *q    = NULL ;

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
    printf(" %e %e\n", x_i[i], y_i[i]) ;
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
    printf("%d %d %d %e %e %e\n",type[i], n1[i], n2[i], E[i], A[i], I[i]) ;
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
    printf("%d %d %e\n",d_n[i], d_d[i], d_v[i]) ;
  }
  fprintf(stderr,"  Have %d supports.\n",n_nfors);


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
    printf("%d %d %e\n",f_n[i], f_d[i], f_v[i]) ;
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
    printf("%d %d %e %e\n",l_e[i], l_d[i], l_v1[i], l_v2[i]) ;
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
  if ((F_val   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 
  if ((u_val   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 

  if ((M   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 
  if ((r   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 
  if ((z   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 
  if ((p   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 
  if ((q   = (float *)malloc(3*n_nodes*sizeof(float))) == NULL) { goto memFree ; } 

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
    for (j=0; j<n_eload; j++)
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
  if ((K_val = (float *)malloc(K_len*sizeof(float))) == NULL) { goto memFree ; } 

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

double md_norm_K()
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

double md_vec_norm(a, len)
double *a;
long len;
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

	normA = md_norm_K();
	normB = md_vec_norm(F_val, n);

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
	    fprintf(stderr,"zero value at [%i,%i]: %e\n",i+1,i+1, M[i]);
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
    fprintf(stderr,"    CG iteration: %i/%i\n",i,n);
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

	  normRes = md_vec_norm(r, n);
	  normX   = md_vec_norm(u_val, n);

    if (normRes  <= ((1e-6)*((normA*normX) + normB)) ) 
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
float l ;
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
  tuh = (E*A)/l ;
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
float s ;
float c ;
{
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

/** Transposed transformation matrix */
void tran_t(s, c)
float s ;
float c ;
{
  T[0][0] = c ;
  T[0][1] = (-s) ;
  T[1][0] = s ;
  T[1][1] = c ;

  T[3][3] = c ;
  T[3][4] = (-s) ;
  T[4][3] = s ;
  T[4][4] = c ;
}

void ke_to_keg(s, c)
float s ;
float c ;
{
  int i,j,k;
  float fval, kval ;

  tran_t(s, c);

  for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      fval += T[j][i] * fe[j] ;

      kval = 0.0 ;
      for (k=0; k<6; k++)
      {
        /* normally: val += T[i][k] * Ke_loc[k][j] ; */
        kval += T[k][i] * ke[k][j] ;
      }
      keg[i][j] = kval ;
    }
  }

  tran(s, c);
  ke_switch();

  for (i=1; i<=6; i++)
  {
    fval = 0.0 ;
    for (j=1; j<=6; j++)
    {
      kval = 0.0 ;
      for (k=1; k<=6; k++)
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

/* computes stiffness matrix of the structure */
void stiff()
{
  int i, j, k, ii, jj ;
  int type ;
  float x1,y1, x2,y2, l, s, c ;


  for (i=0; i<n_elems; i++)
  {
    x1 = x_i[n1[i]-1] ;
    y1 = y_i[n1[i]-1] ;
    x2 = x_i[n2[i]-1] ;
    y2 = y_i[n2[i]-1] ;
    
    l = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) ;
    if (l <= 0.0) {continue;} /* oops, zero length element */
    s = (y2-y1)/l ;
    c = (x2-x1)/l ;
    
    tran_zero();
    stiff_loc(type, E[i], A[i], I[i], l) ;
    ke_to_keg(s, c) ;

    /* TODO: localisation */
    for (k=0; k<6; k++)
    {
      if (k <3) { ii = n1[i]*3 + k - 2 ; }
      else      { ii = n2[i]*3 + k - 5 ; }
  
      for (j=0; j<6; j++)
      {
        if (j <3) { jj = n1[i]*3 + j - 2 ; }
        else      { jj = n2[i]*3 + j - 5 ; }
  
        md_K_add(ii, jj, keg[k][j]) ;
      }
    }
  }
}

/*TODO: forces and boundary conditions! */

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

int main(argc, argv)
int argc ;
char *argv[];
{
  FILE *fw = NULL ;

  if (argc < 2)
  {
    /* standard input */
    fprintf(stderr,"\nInteractive input:\n\n");
    read_data(stdin);
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

  if (alloc_kf() != 0 )
  {
    free_data();
    return(-1);
  }

  stiff(); 

  solve_eqs();

  free_data();
  free_sol_data();

  return(0);
}

/* end of ddfor.c */
