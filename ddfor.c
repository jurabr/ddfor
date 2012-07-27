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

int    n_dofs  = 0 ;

float *x_i = NULL ;
float *y_i = NULL ;

float *E = NULL ;
float *A = NULL ;
float *I = NULL ;
int   *n1 = NULL ;
int   *n2 = NULL ;
int    *type = NULL ;

int   *d_n = NULL ; /* node */
int    *d_d = NULL ; /* direction 1=x 2=y 3=rot */
float *d_v = NULL ; /* size */

int   *f_n = NULL ; /* node */
int    *f_d = NULL ; /* direction 1=fx 2=fy 3=m */
float *f_v = NULL ; /* size */

int   *l_e = NULL ; /* node */
int    *l_d = NULL ; /* direction 1=x 2=y, 3=x global, 4=y global */
float *l_v1 = NULL ; /* size at beginning */
float *l_v2 = NULL ; /* size at end */

/* solution variables: */
float  ke[6][6] ;
float  keg[6][6] ;
float  T[6][6] ;
float  fe[6];
float  feg[6];
float  ueg[6];

float  **K ;
float   *F = NULL ;

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


  n_dofs = n_nodes * 3 ;
  fprintf(stderr,"  Problem size is: %d.\n",n_dofs);

  return(0);
}

/* Allocates space for linear system */
int alloc_kf(dofs, rows)
int dofs;
int rows;
{
  int i,j;

  if ((K = (float **)malloc((dofs+1)*sizeof(float *))) == NULL)
  {
    return(-1);
  }
  else
  {
    K[0] = NULL ;
    for (i=1; i<=dofs; i++)
    {
      if ((K[i] = (float *)malloc((rows+1)*sizeof(float))) == NULL)
      {
        for (j=1; j<i; j++)
        {
          free(K[j]) ; K[j] = NULL ;
        }
        free(K); K = NULL ;
        fprintf(stderr,"No space for K matrix (stopped at %d row)!\n",i);
        return(-2);
      }
      else
      {
        for (j=1; j<=rows; j++)
        {
          K[i][j] = 0.0 ;
        }
      }
    }

    if ((F = (float *)malloc((dofs+1)*sizeof(float))) == NULL)
    {
      for (j=1; j<=dofs; j++)
      {
        free(K[j]) ; K[j] = NULL ;
      }
      free(K); K = NULL ;
      fprintf(stderr,"No space for F matrix!\n");
      return(-2);
    }
  }

  return(0);
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

/* band size computation routine */
int band()
{
  int max = 0 ;
  int i;
	int nind;
	int ni1, ni2;

	for (i=0;i<n_elems; i++)
	{
	   ni1 =n1[i];
	   ni2 =n2[i];
		 if (ni2<ni1)
		 {
		   nind=ni2;
			 ni2=ni1;
			 ni1=nind;
		 }
		 nind=3*(ni2-1)+3-3*(ni1-1);
		 if (max < nind) max=nind;
	}
  fprintf(stderr,"  MAX: %d\n",max);
	if (max == 0)
	{
	  /* if failed */
		fprintf(stderr,"Band computing failed");
	}
	else
	{
    /*
	  max *= 2 ;
		max += 1 ;
    */
	}
	return(max);
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

void stiff()
{
  int i ;
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
  }
}

/** Gauss elimination */
int solve_band(K,F, k_band)
float **K;
float  *F;
int     k_band;
{
  int i,j,k, r,h,g;
  float eps ;

  r = (int)((k_band+1)/2) ;
  eps = 1e-6 ;

  /* L-U decomposition routine: */
  for (k=1; k<=n_dofs; k++)
  {
    if (fabs(K[k][r]) <= eps) {break;}
    K[k][r] = 1.0 / K[k][r] ;
    h = r - 1 ;
    i = k + 1 ;
    if ((h>=1)||(i <= n_dofs))
    {
      K[i][h] = K[i][h] * K[k][r] ;
      j = h + 1 ;
      g = r + 1 ;
      
      if ( (g<=k_band) || (j<=(r+n_dofs-1)) )
      {
        K[i][j] = K[i][j] - (K[i][h]*K[i][g]) ;
        j = j + 1;
        g = g + 1;
      }
      i = i + 1;
      h = h - 1;
    }
  }

  /* forward routine: */
  for (k=1; k<=(n_dofs-1); k++)
  {
    i = k + 1 ;
    j = r - 1 ;

    if ( (j >= 1) || (i <= n_dofs))
    {
      F[i] = F[i] - K[i][j]*F[k] ;
      i = i + 1 ;
      j = j - 1 ;
    }
  }

  /* backward routine: */
  for (k=n_dofs; k>=1; k--)
  {
    i = k + 1 ;
    j = r + 1 ;
    if ((j<=k_band) || (i<=n_dofs))
    {
      F[k] = F[k] - K[k][j]*F[i] ;
      i = i + 1 ;
      j = j + 1 ;
    }
    F[k] = F[k] * K[k][r] ;
  }

  return(0);
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

/** free K,F */
void free_sol_data()
{
  int j ;
  for (j=1; j<=n_dofs; j++)
  {
    free(K[j]) ; K[j] = NULL ;
  }
  free(K); K = NULL ;
  free(F); F = NULL ;
}

int main(argc, argv)
int argc ;
char *argv[];
{
  FILE *fw = NULL ;
  int   n_band = 0 ;

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

  n_band = band() ;
  printf("Matrix band size is %d, total size is %d.\n", n_band,n_nodes*3);

  if (alloc_kf(n_dofs, n_band) != 0 )
  {
    free_data();
    return(-1);
  }

  stiff(); 


  free_data();
  free_sol_data();

  return(0);
}

/* end of ddfor.c */
