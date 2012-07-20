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
float   *u = NULL ;

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
    free(f_n);free(f_d);free(f_v);
    free(l_e);
    fprintf(stderr,"Can't allocate types of loads!\n");
    return(-2);
  }
  if ((l_v1=(float *)malloc(n_eload*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    free(f_n);free(f_d);free(f_v);
    free(l_e);free(l_d);
    fprintf(stderr,"Can't allocate sizes of loads!\n");
    return(-2);
  }
  if ((l_v2=(float *)malloc(n_eload*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    free(d_n);free(d_d);free(d_v);
    free(f_n);free(f_d);free(f_v);
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

/* test transformation matrix to zero */
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
  float E, A, I ;

  for (i=0; i<n_elems; i++)
  {
    /* TODO: read x,y ; */
    /* TODO: read type, E, A, I ; */
    l = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)) ;
    if (l <= 0.0) {continue;} /* ops, zero length element */
    s = (y2-y1)/l ;
    c = (x2-x1)/l ;
    
    tran_zero();
    stiff_loc(type, E, A, I, l) ;
    ke_to_keg(s, c) ;

    /* TODO: localisation */
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
      fprintf(stderr,"FAiled to open input file!\n");
      return(-1);
    }
    else
    {
      read_data(fw);
      fclose(fw);
    }
  }

  stiff(); 

  return(0);
}

/* end of ddfor.c */
