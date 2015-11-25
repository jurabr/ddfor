/*
 * File name: ddfor.c
 * Date:      2012/07/17 20:10
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

int    sol_mode = 0 ; /* 0=statics, 1=stability, 2=modal, 3=aaem */

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
float *rho = NULL ;
int   *n1 = NULL ;
int   *n2 = NULL ;
int   *type = NULL ;
int   *e_g = NULL ; /* element group */

int   *d_n = NULL ; /* node */
int   *d_d = NULL ; /* direction 1=x 2=y 3=rot */
float *d_v = NULL ; /* size */
int   *d_g = NULL ; /* load group */

int   *f_n = NULL ; /* node */
int   *f_d = NULL ; /* direction 1=fx 2=fy 3=m */
float *f_v = NULL ; /* size */
int   *f_g = NULL ; /* load group */

int   *l_e = NULL ; /* node */
int   *l_d = NULL ; /* direction 1=x 2=y, 3=x global, 4=y global */
float *l_v1 = NULL ; /* size at beginning */
float *l_v2 = NULL ; /* size at end */
int   *l_g  = NULL ; /* load group */

/* solution variables: */
double  ke[6][6] ;
double  keg[6][6] ;
double  T[6][6] ;
double  fe[6];
double  feg[6];
double  ueg[6];
double  ue[6];

int    K_len    = 0 ; /* actual size of K and M matrices */
int    K_size   = 0 ; /* number of unknowns              */
double *K_val    = NULL ;
double *F_val    = NULL ;
double *u_val    = NULL ;

double *M_val    = NULL ;
double *Mu_val    = NULL ;
double *Fr_val    = NULL ;
double *uu_val    = NULL ;
double *EE      = NULL ; /* initial value for viscoelastic computations */
double *E1      = NULL ; /* E(t1) for viscoelastic computations */

double *M    = NULL ;
double *r    = NULL ;
double *z    = NULL ;
double *p    = NULL ;
double *q    = NULL ;

int  K_nfree = 0    ; /* number of free rotations */
int *K_fe    = NULL ; /* free rotations elements  */
int *K_fn    = NULL ; /* free rotations nodes     */
int *K_fpos  = NULL ; /* free rotations positions */
int pvec[6]         ; /* localisation field       */

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

  /* allocate data for nodes */
  if ((x_i=(float *)malloc(n_nodes*sizeof(float))) == NULL)
  {
    fprintf(stderr,"Can't allocate X coordinates!\n");
    return(-2);
  }

  if ((y_i=(float *)malloc(n_nodes*sizeof(float))) == NULL)
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
    fscanf(fw,"%e %e", &x_i[i], &y_i[i]) ;
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
  if ((rho=(float *)malloc(n_elems*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);
    fprintf(stderr,"Can't allocate moments of inertia!\n");
    return(-2);
  }

  if ((type=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);
    fprintf(stderr,"Can't allocate types!\n");
    return(-2);
  }
  if ((e_g=(int *)malloc(n_elems*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);
    fprintf(stderr,"Can't allocate types!\n");
    return(-2);
  }


  /* read element data */
#ifdef UI
  fprintf(stderr,"Element data (type node1 node2 E A I dens grp):\n");
#endif
  for (i=0; i<n_elems; i++)
  {
    fscanf(fw,"%d %d %d %e %e %e %e %d",
      &type[i],&n1[i],&n2[i],&E[i],&A[i],&I[i],&rho[i],&e_g[i]) ;
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
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    fprintf(stderr,"No supports!\n");
    return(-1);
  }

  if ((d_n=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    fprintf(stderr,"Can't allocate nodes for supports!\n");
    return(-2);
  }
  if ((d_d=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);
    fprintf(stderr,"Can't allocate types of supports!\n");
    return(-2);
  }
  if ((d_v=(float *)malloc(n_disps*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);
    fprintf(stderr,"Can't allocate sizes of supports!\n");
    return(-2);
  }
  if ((d_g=(int *)malloc(n_disps*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
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
    fscanf(fw,"%d %d %e %d",&d_n[i], &d_d[i], &d_v[i], &d_g[i]) ;
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
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    fprintf(stderr,"Can't allocate nodes for forces!\n");
    return(-2);
  }
  if ((f_d=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    free(f_n);
    fprintf(stderr,"Can't allocate types of forces!\n");
    return(-2);
  }
  if ((f_v=(float *)malloc(n_nfors*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    free(f_n);free(f_d);
    fprintf(stderr,"Can't allocate sizes of forces!\n");
    return(-2);
  }
  if ((f_g=(int *)malloc(n_nfors*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
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
    fscanf(fw,"%d %d %e %d",&f_n[i], &f_d[i], &f_v[i], &f_g[i]) ;
  }
#ifdef UI
  fprintf(stderr,"  Have %d forces in nodes.\n",n_nfors);
#else
  fprintf(stderr,"  Forces in nodes: %d\n",n_nfors);
#endif
  } /* end of n_nfors > 0 */


  /* loads on elements */
#ifdef UI
  fprintf(stderr,"Number of element loads:\n");
#endif
  fscanf(fw,"%d", &n_eload);
  if (n_eload <= 0)
  {
    fprintf(stderr,"No element loads.\n");
  }
  else
  {
  if ((l_e=(int *)malloc(n_eload*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);
    free(f_n);free(f_d);free(f_v);
    fprintf(stderr,"Can't allocate elements for loads!\n");
    return(-2);
  }
  if ((l_d=(int *)malloc(n_eload*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);free(f_g);}
    free(l_e);
    fprintf(stderr,"Can't allocate types of loads!\n");
    return(-2);
  }
  if ((l_v1=(float *)malloc(n_eload*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);free(f_g);}
    free(l_e);free(l_d);
    fprintf(stderr,"Can't allocate sizes of loads!\n");
    return(-2);
  }
  if ((l_v2=(float *)malloc(n_eload*sizeof(float))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);free(f_g);}
    free(l_e);free(l_d);free(l_v1);
    fprintf(stderr,"Can't allocate sizes of loads!\n");
    return(-2);
  }
  if ((l_g=(int *)malloc(n_eload*sizeof(int))) == NULL)
  {
    free(x_i); free(y_i); free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
    free(d_n);free(d_d);free(d_v);free(d_g);
    if (n_nfors > 0){ free(f_n);free(f_d);free(f_v);free(f_g);}
    free(l_e);free(l_d);free(l_v1);free(l_v2);
    fprintf(stderr,"Can't allocate sizes of loads!\n");
    return(-2);
  }

  
  /* read element loads data */
#ifdef UI
  fprintf(stderr,"Element loads (element direction startsize endsize grp):\n");
#endif
  for (i=0; i<n_eload; i++)
  {
    fscanf(fw,"%d %d %e %e %d",&l_e[i], &l_d[i], &l_v1[i], &l_v2[i],&l_g[i]) ;
  }
#ifdef UI
  fprintf(stderr,"  Have %d element loads.\n",n_eload);
#else
  fprintf(stderr,"  Element loads:   %d\n",n_eload);
#endif
  } /* end of n_eloads > 0*/

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
    fprintf(fw,"%d %d %d %e %e %e %e %d\n",
      type[i],n1[i],n2[i],E[i],A[i],I[i],rho[i],e_g[i]) ;
  }

  fprintf(fw,"%d\n", n_disps);
  for (i=0; i<n_disps; i++)
  { fprintf(fw,"%d %d %e %d\n",d_n[i], d_d[i], d_v[i], d_g[i]) ; }

  fprintf(fw,"%d\n", n_nfors);
  for (i=0; i<n_nfors; i++)
  { fprintf(fw,"%d %d %e %d\n",f_n[i], f_d[i], f_v[i], f_g[i]) ; }

  fprintf(fw,"%d\n", n_eload);
  for (i=0; i<n_eload; i++)
  { fprintf(fw,"%d %d %e %e %d\n",l_e[i], l_d[i], l_v1[i], l_v2[i], l_g[i]) ; }

  return(0);
}
#endif

/** free K,F */
void free_sol_data()
{
  if (F_val != NULL )free(F_val);
  if (u_val != NULL )free(u_val);
  if (K_val != NULL )free(K_val);

  if (sol_mode > 0)
  {
    if (M_val !=NULL)free(M_val);
    if (Mu_val !=NULL)free(Mu_val);
    if (Fr_val !=NULL)free(Fr_val);
    if (uu_val !=NULL)free(uu_val);
    if (sol_mode == 3)
    {
      if (EE !=NULL)free(EE);
      if (E1 !=NULL)free(E1);
    }
  }
}

/* Compute free rotations */
void comp_frot()
{
  int i ;
  int sum ;

  sum = 0 ;

  for (i=0; i<n_elems; i++)
  {
    switch (type[i])
    { 
      case 3: sum = sum + 2 ; break ;
      case 2: sum = sum + 1 ; break ;
      case 1: sum = sum + 1 ; break ;
      case 0: 
      default: break ;
    }
  }

  if (sum > 0)
  {
    if ((K_fe = (int *)malloc(sum*sizeof(int)))   == NULL) { goto memFree;}
    if ((K_fn = (int *)malloc(sum*sizeof(int)))   == NULL) { goto memFree;}
    if ((K_fpos = (int *)malloc(sum*sizeof(int)))   == NULL) { goto memFree;}
    K_nfree = sum ;
  }
  else
  {
    return ;
  }

  sum = 0 ;
  for (i=0; i<n_elems; i++)
  {
    switch (type[i])
    { 
      case 3: 
              K_fe[sum]   = i ;
              K_fn[sum]   = 1 ;
              K_fpos[sum] = (3*n_nodes) + sum + 1 ;
              sum = sum + 1 ;

              K_fe[sum]   = i ;
              K_fn[sum]   = 2 ;
              K_fpos[sum] = (3*n_nodes) + sum + 1 ;
              sum = sum + 1 ; 
              break ; /* o--o */
      case 2:
              K_fe[sum]   = i ;
              K_fn[sum]   = 2 ;
              K_fpos[sum] = (3*n_nodes) + sum + 1 ;
              sum = sum + 1 ; 
              break ; /* |--o */
      case 1: 
              K_fe[sum]   = i ;
              K_fn[sum]   = 1 ;
              K_fpos[sum] = (3*n_nodes) + sum + 1 ;
              sum = sum + 1 ;
              break ; /* o--| */
      case 0: 
      default: break ;
    }
  }
  return ;
memFree:
  free(K_fe); K_fe = NULL ;
  free(K_fn); K_fn = NULL ;
  free(K_fpos); K_fpos = NULL ;
  K_nfree = 0 ;
}

/** fills localisation vector for free roation on element */
void e_frotv(epos)
  int epos ;
{
  int i ;

  for (i=0; i<6; i++) /* normal routine */
  {
    if (i <3) { pvec[i] = n1[epos]*3 + i - 2 ; }
    else      { pvec[i] = n2[epos]*3 + i - 5 ; }
  }

  for (i=0; i<K_nfree; i++)
  {
    if (epos == K_fe[i])
    {
      switch (type[epos])
      {
        case 3: /* o--o */
          pvec[2] = K_fpos[i] ; pvec[5] = K_fpos[i+1] ; return ; break;
          break;
        case 2: /* |--o */
          pvec[5] = K_fpos[i] ; return ; break;
          break;
        case 1:  /* o--| */
          pvec[2] = K_fpos[i] ; return ; break;
        default: /* nothing to do */
          break;
      }
    }
  }
}

/* Allocates space for linear system */
int alloc_kf()
{
  int i ;
  int n = 0 ;

  comp_frot() ; /* computation of free rotations */

  n  = 3*n_nodes + K_nfree ;
  K_len = (n*n - n)/2 + n ;
  K_size = n ;


  if ((F_val = (double *)malloc(K_size*sizeof(double))) == NULL) {goto memFree;} 
  if ((u_val = (double *)malloc(K_size*sizeof(double))) == NULL) {goto memFree;} 

  for (i=0; i<K_size; i++) 
  { 
    F_val[i]   = 0.0 ; 
    u_val[i]   = 0.0 ; 
  } 

  if ((K_val = (double *)malloc(K_len*sizeof(double))) == NULL) { goto memFree ; } 
  for (i=0; i<K_len; i++) { K_val[i]  = 0.0 ; }

  if (sol_mode > 0)
  {
    if ((M_val = (double *)malloc(K_len*sizeof(double))) == NULL) { goto memFree ; } 
    if ((Mu_val = (double *)malloc(K_size*sizeof(double))) == NULL) { goto memFree ; } 
    if ((Fr_val = (double *)malloc(K_size*sizeof(double))) == NULL) { goto memFree ; } 
    if ((uu_val = (double *)malloc(K_size*sizeof(double))) == NULL) { goto memFree ; } 
    for (i=0; i<K_len; i++) { M_val[i]  = 0.0 ; }
    for (i=0; i<K_size; i++)
    {
      Mu_val[i] = 0.0 ;
      Fr_val[i] = 0.0 ;
      uu_val[i] = 0.0 ;
    }
    if (sol_mode == 3)
    {
      if ((EE = (double *)malloc(n_elems*sizeof(double))) == NULL) { goto memFree ; } 
      if ((E1 = (double *)malloc(n_elems*sizeof(double))) == NULL) { goto memFree ; } 
    }
  }

  return(0);
memFree:
  fprintf(stderr,"Not enough memory!");
  free_sol_data();
  return(-1);
}


/* Equation solver (symmetric Gauss): */
int solve_eqs(A, x, b, n)
double *A;
double *x;
double *b;
int n;
{
  int    i,j,k,l,m;
  double c, s;

  for (i=0; i<n; i++) {x[i] = b[i] ;}

  for (k=1; k<n; k++)
  {
    for (i=k+1; i<=n; i++)
    {
      l = k + (i*i-i)/2-1 ;
      m = k + (k*k-k)/2-1 ;
      if (fabs(A[m]) < 1e-6) {return(-1);}
      c = A[l]/A[m] ;
      for (j=i; j<=n; j++)
      {
        l = i + (j*j-j)/2-1 ;
        m = k + (j*j-j)/2-1 ;
        A[l] -= c*A[m] ;
      }
      x[i-1] = x[i-1] - c*x[k-1] ;
    }
  }
  l = n + (n*n-n)/2-1 ;
  x[n-1] = x[n-1] / A[l] ;

  for (i=n-1 ; i >=1; i--)
  {
    s = 0.0 ;
    for (j=i+1; j<=n; j++)
    {
      l = i + (j*j-j)/2-1 ;
      s += A[l] * x[j-1] ;
    }
    l = i + (i*i-i)/2-1 ;
    x[i-1] = (x[i-1]-s) / A[l] ;
  }

  return(0);
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
}

/** Local geometric stiffness matrix */
void geom_loc(type, N, L)
int type ;
float N ;
double L ;
{
  int    i, j ;

  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      ke[i][j] = 0.0 ;
      keg[i][j] = 0.0 ;
    }
  }

  ke[1][1] = 6.0 / 5.0 ;
  ke[1][2] = L / 10.0 ;
  ke[1][4] = -6.0 / 5.0 ;
  ke[1][5] = L / 10.0 ;

  ke[2][1] = L / 10.0 ;
  ke[2][2] = (2*L*L) / 15.0 ;
  ke[2][4] = -L / 10.0 ;
  ke[2][5] = -(L*L) / 30.0 ;

  ke[4][1] = -6.0 / 5.0 ;
  ke[4][2] = -L / 10.0 ;
  ke[4][4] = 6.0 / 5.0 ;
  ke[4][5] = -L / 10.0 ;

  ke[5][1] = L / 10.0 ;
  ke[5][2] = -(L*L) / 30.0 ;
  ke[5][4] = -L / 10.0 ;
  ke[5][5] = (2*L*L) / 15.0 ;

  ke[0][0] = 1.0 ;
  ke[0][3] = -1.0 ;
  ke[3][0] = -1.0 ;
  ke[3][3] = 1.0 ;

  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++) { ke[i][j] = (N/L) * ke[i][j] ; }
  }
}

/** Local mass matrix */
void mass_loc(type, dens, Ax, Lx)
int type ;
float dens ;
float Ax ;
double Lx ;
{
  int i, j ;
  float mass ;

  for (i=0; i<6; i++)
  {
    for (j=0; j<6; j++)
    {
      ke[i][j] = 0.0 ;
      keg[i][j] = 0.0 ;
    }
  }
	mass = dens * Ax * Lx / 420.0 ;
  ke[0][0] = mass*140.0 ;
  ke[0][3] = mass*70.0 ;

	ke[1][1] = mass*156.0 ;
  ke[1][2] = 22.0*Lx*mass ;
  ke[1][4] = 54.0*mass ;
  ke[1][5] = (-13.0)*Lx*mass ;

  ke[2][2] = mass * 4.0*Lx*Lx ;
  ke[2][4] = mass * 13.0*Lx ;
  ke[2][5] = mass * (-3.0)*Lx*Lx ;

  ke[3][3] = 140.0*mass ;

  ke[4][4] = mass * 156.0 ;
  ke[4][5] = mass* (-22.0)*Lx ;

  ke[5][5] = mass * 4.0*Lx*Lx ;

	for (i=0; i<6; i++)
	{
		for (j=i; j<6; j++) { ke[j][i] =  ke[i][j] ; }
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

/* fills "ke" with content of "keg" */
void ke_switch()
{
  int i, j;
  for (i=0; i<6; i++) { for (j=0; j<6; j++) { ke[i][j] = keg[i][j] ; } }
}

/* tranforms ke and fe to keg and feg (set fev=0 to ignore fe->feg)*/
void ke_to_keg(s, c, fev)
double s ;
double c ;
int fev ;
{
  int i,j,k;
  float fval, kval ;

  tran(s, c);

  for (i=0; i<6; i++)
  {
    fval = 0.0 ;

    for (j=0; j<6; j++)
    {
      if (fev != 0) fval += T[j][i] * fe[j] ;

      kval = 0.0 ;
      for (k=0; k<6; k++)
      {
        kval += T[k][i] * ke[k][j] ;
      }
      keg[i][j] = kval ;
    }
    if (fev != 0)feg[i] = fval ;
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
  if (col < row) return  ;
  i = row + (col*col-col)/2-1 ;
  K_val[i] += val ;
}

/* puts data to the right place in M */
void md_M_add(row, col, val)
int row;
int col;
float val;
{
  int i ;
  if (col < row) return  ;
  i = row + (col*col-col)/2-1 ;
  M_val[i] += val ; 

}


/* loads on elements */
void one_eload(epos, na, nb, va, vb, L)
int epos;
double na;
double nb;
double va;
double vb;
double L;
{
  fe[0]+=(-(2.0*na+1.0*nb)*L)/6.0 ;
  fe[1]+=(-(7.0*va+3.0*vb)*L)/20.0 ;
  fe[2]+=((3.0*va+2.0*vb)*L*L)/60.0 ;
  fe[3]+=(-(1.0*na+2.0*nb)*L)/6.0 ;
  fe[4]+=(-(3.0*va+7.0*vb)*L)/20.0 ;
  fe[5]+=(-(2.0*va+3.0*vb)*L*L)/60.0 ;
}

/* computes stiffness matrix of the structure */
void stiff(eg, lc)
int eg;
int lc;
{
  int i, j, k, ii, jj, m ;
  float x1,y1, x2,y2, l, s, c ;

  for (i=0; i<n_elems; i++)
  {
    if ((e_g[i] > 0)&&(e_g[i] != eg)) continue;
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
    for (m=0; m<n_eload; m++)
    {
      if ((l_g[i] != lc)) continue; /* only lc-th step is considered */
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
    ke_to_keg(s, c,1) ;

    /* localisation */
    e_frotv(i, pvec);
    for (k=0; k<6; k++)
    {
      ii = pvec[k] ;
      F_val[ii-1] += feg[k] ;

      for (j=0; j<6; j++)
      {
        jj = pvec[j] ;
        md_K_add(ii, jj, keg[k][j]) ;
      }
    }

    /* mass matrix */
    if (sol_mode == 2) /* modal */
    {
      tran_zero();
      mass_loc(type, rho[i], A[i], (double)l);
      ke_to_keg(s, c,0) ;
      for (k=0; k<6; k++)
      {
        ii = pvec[k] ;
        for (j=0; j<6; j++)
        {
          jj = pvec[j] ;
          md_M_add(ii, jj, keg[k][j]) ;
        }
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

  row = (node-1)*3 + dir  ;  /* note: -1 ? */
  n = K_size ;

  for (i=1; i<=n; i++)
  {
    if (i > row) { j = row + (i*i-i)/2-1 ; }
    else         { j = i + (row*row-row)/2-1 ; }
    F_val[i-1] += K_val[j]*val ;
  }

  for (i=1; i<=n; i++)
  {
    if (i > row) { j = row + (i*i-i)/2-1 ; }
    else         { j = i + (row*row-row)/2-1 ; }
    K_val[j] = 0.0 ;
  }

  K_val[row+(row*row-row)/2-1] = 1.0 ;

  u_val[row-1] = val ;
  F_val[row-1] = val ; /* it will destroy any force in this place */
}

/* adds supports to M (zero values only) */
void add_one_disp_M(node, dir)
int   node;
int   dir;
{
  int i,j,n ;
  int row ;

  row = (node-1)*3 + dir  ;  /* note: -1 ? */
  n = K_size ;

  for (i=1; i<=n; i++)
  {
    if (i > row) { j = row + (i*i-i)/2-1 ; }
    else         { j = i + (row*row-row)/2-1 ; }
    M_val[j] = 0.0 ;
  }

  M_val[row+(row*row-row)/2-1] = 1.0 ;
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

  if (sol_mode > 0) 
  {
    for (i=0; i<n_disps; i++) 
    { 
      if ((d_g[i] > 0)||(d_g[i] != lc)) continue;
      add_one_disp_M(d_n[i], d_d[i]); 
    } 
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
    free(n1); free(n2);free(E);free(A);free(I);free(rho);free(type);free(e_g);
  }
  if (n_disps > 0)
  {
    free(d_n);free(d_d);free(d_v);free(d_g);
  }
  if (n_nfors > 0)
  {
    free(f_n);free(f_d);free(f_v);free(f_g);
  }
  if (n_eload > 0)
  {
    free(l_e);free(l_d);free(l_v1);free(l_v2);free(l_g);
  }
  if (K_nfree > 0)
  {
    free(K_fe); K_fe = NULL ;
    free(K_fn); K_fn = NULL ;
    free(K_fpos); K_fpos = NULL ;
    K_nfree = 0 ;
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

/** Find all element load in given direction */
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

/* Element deformation computation:  */
double in_def(e_type, L, q1, q2, x)
int e_type;
double L;
double q1;
double q2;
double x;
{
  double y = 0.0 ;
  double xx;
  double qa, qb ;

  if (q1 == q2) { qa = q1 ; qb = 0.0 ; }
  else          { qa = q1 ; qb = q2 - q1 ; }

  switch(e_type)
  {
    case 0: /* |-| */
            y  = (qa*x*x*x*x)/24 - (qa*L*x*x*x)/12 + (qa*L*L*x*x)/24 ;
            y += ((qb*x*x)/(120*L))*(x*x*x - 3*L*L*x + 2*L*L*L) ;
            break;
    case 1: /* o-| */ 
            y  = (qa*x*x*x*x)/2 - (3*qa*L*x*x*x)/48 + (qa*L*L*L*x)/8 ;
            y += (qb*x*x*x*x)/24 - (3*qb*L*x*x*x)/48 + (qb*L*L*L*x)/48 ;
            break;
    case 2: /* |-o */ 
            y  = (qa*x*x*x*x)/2 - (3*qa*L*x*x*x)/48 + (qa*L*L*L*x)/8 ;
            xx = L-x ;
            y += (qb*xx*xx*xx*xx)/24 - (3*qb*L*xx*xx*xx)/48 + (qb*L*L*L*xx)/48 ;
            break;
    case 3: /* o-o */ 
            y  = (qa*x*x*x*x)/24 - (qa*L*x*x)/12 + (qa*L*L*L*x)/24 ;
            y += (qb*x*x*x*x*x)/(120*L) - (qb*L*x*x*x)/36 + (7/360)*qb*L*L*L*x ;
            break;
  }
  return((-1.0)*y);
}

/** Compute internal force (N, V, M) for given point of beam */
double in_force(ftyp, epos, div, ppos)
int ftyp; /* force type: 1=N, 2=V, 3=M, 4=w */
int epos; /* element position */
int div;  /* number of divisions */
int ppos; /* number of computed point (0...div)*/
{
  double x1,x2,y1,y2 ;
  double Na,Ma,Mb, L, lenx, lenxx, na, nb, no, nt ;
  double Xo = 0.0 ;

  Na = fe[0];
  Ma = fe[2];
  Mb = fe[5];

  x1 = x_i[n1[epos]-1] ;
  y1 = y_i[n1[epos]-1] ;
  x2 = x_i[n2[epos]-1] ;
  y2 = y_i[n2[epos]-1] ;
 

  L = sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) ) ;
  lenx  = L*((double)((double)ppos/(double)(div))) ;
  lenxx = L - lenx ;

  switch (ftyp)
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
    case 4 : get_eloads(epos, 2, &na, &nb);
             return(in_def(type[epos], L, na, nb, lenx));
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
  e_frotv(epos);
  for (i=0; i<3; i++)
  {
    ueg[i]   = u_val[pvec[i]-1];
    ueg[i+3] = u_val[pvec[i+3]-1];
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
  for (m=0; m<n_eload; m++)
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

  fprintf(fw,"\nDDFOR FINAL REPORT\n");

  fprintf(fw,"\nNodes:\n");
  fprintf(fw," Num     X            Y:\n");
  for (i=0; i<n_nodes; i++)
  {
    fprintf(fw," %3d %e %e\n",i+1, x_i[i], y_i[i]);
  }

  fprintf(fw,"\nElements:\n");
  fprintf(fw," Num Node1 Node2     E             A           I:            Grp:\n");
  for (i=0; i<n_elems; i++)
  {
    fprintf(fw," %3d   %3d   %3d %e %e %e ",i+1,n1[i],n2[i],E[i],A[i],I[i]);
    switch(type[i])
    {
      case 0: fprintf(fw,"|--|");break;
      case 1: fprintf(fw,"o--|");break;
      case 2: fprintf(fw,"|--o");break;
      case 3: fprintf(fw,"o--o");break;
      default: fprintf(fw,"unknown!");break;
    }
    fprintf(fw," %2d\n",e_g[i]);
  }

  fprintf(fw,"\nSupports in nodes:\n");
  fprintf(fw," Num Node Dir Size:        LC:\n");
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
    fprintf(fw," %e %3d\n",d_v[i],d_g[i]);
  }

  fprintf(fw,"\nLoads in nodes:\n");
  fprintf(fw," Num Node Dir Size:           LC:\n");
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
    fprintf(fw," %e %3d\n",f_v[i],f_g[i]);
  }

  fprintf(fw,"\nElement loads:\n");
  fprintf(fw," Num Node Dir Sizes (start..end):         LC:\n");
  for (i=0; i<n_eload; i++)
  {
    fprintf(fw," %3d  %3d ",i+1, l_e[i]);
    switch(l_d[i])
    {
      case 1: fprintf(fw,">>>");break;
      case 2: fprintf(fw,"vvv");break;
      default: fprintf(fw,"unknown!");break;
    }
    fprintf(fw," %e..%e %3d\n",l_v1[i],l_v2[i],l_g[i]);
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

/* geometric stiffness matrix */
void geom_stiff(eg)
int eg;
{
#ifdef LARGE
  int i, j, k, ii, jj, m ;
  float x1,y1, x2,y2, l, s, c, fval ;
  double N ;

  if (sol_mode != 1) return;


  for (i=0; i<n_elems; i++)
  {
    if ((e_g[i] > 0)||(e_g[i] != eg)) continue;
    e_frotv(i);
    for (j=0; j<3; j++)
    {
      ueg[j]   = u_val[pvec[j]-1];
      ueg[j+3] = u_val[pvec[j+3]-1];
    }

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
    u_to_ue(s,c);
    stiff_loc(type[i], E[i], A[i], I[i], (double)l) ;
    /* compute forces kn nodes: */
    for (k=0; k<6; k++)
    {
      fval = 0.0 ;
      for (j=0; j<6; j++) { fval += ke[k][j] * ue[j] ; }
      fe[k] -= fval ;
    }

    N = -0.5*(fe[0] + fe[4]) ; 
    geom_loc(type, (double)N, (double)l);
    ke_to_keg(s, c,0) ;
    for (k=0; k<6; k++)
    {
      ii = pvec[k] ;
      for (j=0; j<6; j++)
      {
        jj = pvec[j] ;
        md_M_add(ii, jj, keg[k][j]) ;
      }
    }
  }
#endif
}

/** Compute internal force (N, V, M) for given point of beam
 *  This routine is meant for plotting tools */
void in_gfx(type, epos, div, ppos, mult, vx, vy)
int type; /* force type: 1=N, 2=V, 3=M, 4=w */
int epos; /* element position */
int div;  /* number of divisions */
int ppos; /* number of computed point (0...div)*/
double mult;
double *vx;
double *vy;
{
  double x1,x2,y1,y2,c,s ;
  double L, lenx ;
  double val = 0.0 ;

  x1 = x_i[n1[epos]-1] ;
  y1 = y_i[n1[epos]-1] ;
  x2 = x_i[n2[epos]-1] ;
  y2 = y_i[n2[epos]-1] ;
 

  L = sqrt( (y2-y1)*(y2-y1) + (x2-x1)*(x2-x1) ) ;
  lenx  = L*((double)((double)ppos/(double)(div))) ;

  s = (y2-y1)/L ;
  c = (x2-x1)/L ;

  *vx = x1 + lenx*c ;
  *vy = y1 + lenx*s ;
  val = 0 ;

  switch (type)
  {
    case 0 : /* already computed */ break;
    case 1 : 
    case 2 :
    case 4 : val = in_force(type, epos, div, ppos);
             break ;
    case 3 : val = (-1.0)*in_force(type, epos, div, ppos);
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
#ifdef LARGE
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
  int j ;
  double pos, val, val0;

  *nmax = 0.0 ;
  *vmax = 0.0 ;
  *mmax = 0.0 ;

  for (j=0; j<=div; j++)
  {
    res_loc(epos);
    in_gfx(0, epos, div, j, 1.0, &pos,  &val0);

    in_gfx(1, epos, div, j, 1.0, &val, &val0);
    if (fabs(val) > fabs(*nmax)) {*nmax = val ; *npos = pos;}

    in_gfx(2, epos, div, j, 1.0, &val, &val0);
    if (fabs(val) > fabs(*vmax)) {*vmax = val ; *vpos = pos;}

    in_gfx(3, epos, div, j, 1.0, &val, &val0);
    if (fabs(val) > fabs(*mmax)) {*mmax = val ; *mpos = pos;}
  }

  return(0);
}
#endif

/** Plot data (model geometry) to text console */
#ifndef NO_PSEUDO_GFX
void pseudo_geom(fw, mode)
FILE *fw;
int   mode;
{
  int size_x = 78 ; /* x console size (width)  */
  int size_y = 24 ; /*  y console size (height) */
  double  mult_x, mult_y ;
  double max_x, min_x, max_y, min_y ;
  double x_e, y_e, dx, dy ;
  char   **fld = NULL;
  char   symbol;
  char   str[5];
  int i, j;
  int ii, jj;
  int is, ilen ;

#ifdef POFO
  size_x = 39 ; size_y = 8 ; /* atari portfolio screen */
#endif

  /* allocate image buffer: */
  if ((fld = (char **)malloc (size_y*sizeof(char *))) == NULL) return;
  for (i=0; i<size_y; i++)
  {
    if ((fld[i] = (char *)malloc (size_x*sizeof(char))) == NULL) return;
    for (j=0; j<size_x; j++) { fld[i][j] = ' '; }
  }

  /* find max, min: */
  max_x = x_i[0] ;
  min_x = x_i[0] ;
  max_y = y_i[0] ;
  min_y = y_i[0] ;

  for (i=1; i<n_nodes; i++)
  {
    if (max_x < x_i[i]) {max_x = x_i[i];}
    if (min_x > x_i[i]) {min_x = x_i[i];}

    if (max_y < y_i[i]) {max_y = y_i[i];}
    if (min_y > y_i[i]) {min_y = y_i[i];}
  }

  /* compute multipliers for translation to screen size: */
  if (fabs(max_x-min_x) > 0)
  {mult_x = (double)(size_x-3) /  (max_x - min_x) ;}
  else { mult_x = 1 ; }

  if (fabs(max_y-min_y) > 0)
  { mult_y = (double)(size_y-3) /  (max_y - min_y) ; }
  else { mult_y = 1 ; }

  /* plot elements with numbers */
  for (i=0; i<n_elems; i++)
  {
    /* dx and dy computation */
    x_e =  (x_i[n2[i]-1] - x_i[n1[i]-1]) ;
    y_e =  (y_i[n2[i]-1] - y_i[n1[i]-1]) ;
    symbol = '*';
    if (fabs(x_e) < 1e-6) symbol = '|' ;
    if (fabs(y_e) < 1e-6) symbol = '-' ;
    if (symbol == '*')
    {
      if ((y_e*x_e) > 0.0) symbol = '/' ;
      else                 symbol ='\\';
    }
    /* gfx lenght of element */
    ii =  (int)(fabs(x_e) * mult_x) + 1  ;
    jj =  (int)(fabs(y_e) * mult_y) + 1 ;
    ilen = (int)(sqrt(pow((float)ii,2) + pow((float)jj,2))) ;
    dx = x_e ; dy = y_e ;
    for (is=1; is<ilen; is++)
    {
      x_e = ((float)is/(float)ilen)*dx + x_i[n1[i]-1] ;
      y_e = ((float)is/(float)ilen)*dy + y_i[n1[i]-1] ;
      ii =  (int)((x_e-min_x) * mult_x) + 1  ;
      jj = size_y - ( (int)((y_e-min_y) * mult_y) + 1 ) ;
      fld[jj][ii] = symbol ;
      
      /* hinges: */
      if ((is==2)&&((type[i]==1)||(type[i]==3))) fld[jj][ii] = 'o' ;
      if ((is==(ilen-2))&&((type[i]==2)||(type[i]==3))) fld[jj][ii] = 'o' ;
    }

    if (mode == -1) /* show element numbers */
    {
      /* element numbers code */
      x_e = 0.5 * (x_i[n2[i]-1] + x_i[n1[i]-1]) ;
      y_e = 0.5 * (y_i[n2[i]-1] + y_i[n1[i]-1]) ;
      ii =  (int)((x_e-min_x) * mult_x) + 1  ;
      jj = size_y - ( (int)((y_e-min_y) * mult_y) + 1 ) ;
  
      for (is=0; is<5; is++) { str[is] = '\0'; }
      sprintf(str,"%d",i+1);
      fld[jj][ii] = str[0] ;
      if ((i+1) > 9) fld[jj][ii+1] = str[1] ;
    }
  }

  /* plot "+" for nodes */
  for (i=0; i<n_nodes; i++)
  {
    ii =  (int)((x_i[i]-min_x) * mult_x) + 1  ;
    jj = size_y - ( (int)((y_i[i]-min_y) * mult_y) + 1 ) ;
    if (mode != 0) { fld[jj][ii] = '+' ; }
    else 
    {
      for (is=0; is<5; is++) { str[is] = '\0'; }
      sprintf(str,"%d",i+1);
      fld[jj][ii] = str[0] ;
      if ((i+1) > 9) fld[jj][ii+1] = str[1] ;
    }
  }

  /* plot symbols for supports: "A"=uy, "<"=ux, "X"=rot */
  for (j=0; j<n_disps; j++)
  {
    i = d_n[j]-1;
    ii =  (int)((x_i[i]-min_x) * mult_x) + 1  ;
    jj = size_y - ( (int)((y_i[i]-min_y) * mult_y) + 1 ) ;
    if (mode != 0) 
    { 
       if (d_d[j]==1) fld[jj][ii+1] = '<' ; 
       if (d_d[j]==2) fld[jj][ii] = 'A' ;
       if (d_d[j]==3) 
       {
          if (fld[jj][ii] == 'A' )
          {
             fld[jj-1][ii] = 'A' ;
             fld[jj][ii] = 'X'  ;
           }
           else  fld[jj][ii] = 'X' ;
       }
    }
  } 

  /* plot symbols for forces */
  for (j=0; j<n_nfors; j++)
  {
    i = f_n[j]-1;
    ii =  (int)((x_i[i]-min_x) * mult_x) + 1  ;
    jj = size_y - ( (int)((y_i[i]-min_y) * mult_y) + 1 ) ;
    if (mode != 0) 
    {
      switch (f_d[j]) /* TODO: sizes, orientation */
      {
        case 1: fld[jj][ii] = '=' ;
                fld[jj][ii+1] = '>' ;
                break ;
        case 2: fld[jj-1][ii] = '|'  ;
                fld[jj][ii] = 'V' ;
                break ;
        case 3: fld[jj][ii] = 'M'  ;
                break ;
      }
    }
  }
 

  /* plot data: */  
  fprintf(fw,"\n");
  for (i=0; i<size_y; i++)
  {
    for (j=0; j<size_x; j++)
    {
      fprintf(fw,"%c",fld[i][j]);
    }
    fprintf(fw,"\n");
  }

  /* free data: */
  for (i=0; i<size_y; i++)
  {
    free(fld[i]) ; fld[i] = NULL ;
  }
  free(fld) ; fld = NULL ;
}
#endif

/** Computation of eigenvalues */
int inv_iter(num_res)
  int num_res;
{
  int rv = -1 ;
#if 0
  int i, j, k, jj ;
  double om_top ;
  double om_bot ;
  double omega = 0.0 ;
  double omega0 = 0.0 ;
  double c, mval ;
  double max_iter = 10 ;
  double **eig_val = NULL ;

  if ((eig_val = (double **)malloc(num_res*sizeof(double *))) == NULL)
  {
    fprintf(stderr,"Out of memory for eigen data!\n"); return(rv);
  }
  for (i=0; i<num_res; i++)
  {
    if ((eig_val[i] = (double *)malloc((K_size*sizeof(double)))) == NULL)
    {
      fprintf(stderr,"Out of memory for eigen data!\n"); return(rv);
    }
    for (j=0; j<K_size; j++) { eig_val[i][j] = 0.0 ; }
  }

  omega0  = 0.0 ;
  omega   = 0.0 ;

  for (k=0; k<K_size; k++) /* Mu = M*u ... initial z1 */
  {
    mval = 0.0 ;
    for (j=0; j<K_sizes[k]; j++)
    {
      if  (K_cols[K_from[k]+j] < 0) {break;}
      mval += M_val[K_from[k]+j] * u_val[K_cols[K_from[k]+j]];
    }
    Mu_val[k] = mval ;
  }

  for (j=1; j<=num_res; j++) /* main loop*/
  {
    /* initial approx. */
    if (j > 1) {for(jj=0;jj<K_size;jj++){u_val[jj] *= rand()*1000;}}

    for (i=0; i<max_iter; i++)
    {
      if (j > 1)
      {
        /* Gram-Schmidt: */
        for(jj=0;jj<K_size;jj++)
        {
          F_val[jj] = 0.0;
          Fr_val[jj] = 0.0;
          Mu_val[jj] = 0.0;
        }
        for (k=0; k<K_size; k++) /* Mu = M*u ... initial z1 */
        {
          mval = 0.0 ;
          for (jj=0; jj<K_sizes[k]; jj++)
          {
            if  (K_cols[K_from[k]+jj] < 0) {break;}
            mval += M_val[K_from[k]+jj] * u_val[K_cols[K_from[k]+jj]];
          }
          Mu_val[k] = mval ;
        }
        for (jj=0; jj<(j-1); jj++)
        {
          c = 0.0 ;
          for (k=0; k<K_size; k++) { c += Mu_val[i]*eig_val[jj][k] ; }
          for (k=0; k<K_size; k++) { Fr_val[k] += (c*eig_val[jj][k]) ; }
        }
        for (k=0; k<K_size; k++)
        {
          F_val[k] = u_val[k] - Fr_val[k] ;
          Mu_val[k] = 0.0 ;
        }
        for (k=0; k<K_size; k++) /* Mu = M*F ... z1 for j>=2 */
        {
          mval = 0.0 ;
          for (jj=0; jj<K_sizes[k]; jj++)
          {
            if  (K_cols[K_from[k]+jj] < 0) {break;}
            mval += M_val[K_from[k]+jj] * F_val[K_cols[K_from[k]+jj]];
          }
          Mu_val[k] = mval ;
        }
      }
      for (k=0; k<K_size; k++) /* switch data for solve_eqs: Mu -> F */
      {
        mval = F_val[k] ;
        F_val[k]  = Mu_val[k] ;
        Mu_val[k] = mval ;
      }
      solve_eqs(K_val, u_val, F_val, K_size); /* equation solver*/
      for (k=0; k<K_size; k++) /* switch data for solve_eqs F -> Mu */
      {
        mval = F_val[k] ;
        F_val[k]  = Mu_val[k] ;
        Mu_val[k] = mval ;
      }

      for (k=0; k<K_size; k++) /* Mu = uu (eig_x) */
      {
        mval = 0.0 ;
        for (jj=0; jj<K_sizes[k]; jj++)
        {
          if  (K_cols[K_from[k]+jj] < 0) {break;}
          mval += M_val[K_from[k]+jj] * u_val[K_cols[K_from[k]+jj]];
        }
        uu_val[k] = mval ;
      }
      for (k=0; k<K_size; k++) { om_top += u_val[k]*Mu_val[k] ; }
      for (k=0; k<K_size; k++) { om_bot += u_val[k]*uu_val[k] ; }
      if (fabs(om_bot) <(1e-8))
      {
        fprintf(stderr,"\nInverse iteration failed!\n");
        rv = -1;
        goto memFree;
      }
      omega = om_top / om_bot ; /* eigenvalue */
      printf("OMEGA[%i]: %e\n",i+1,omega);

      for(jj=0;jj<K_size;jj++) {Mu_val[jj]=(uu_val[jj]/sqrt(om_bot)); } /* z(k+1) */

      /* eigenvector (u): */
      for(jj=0;jj<K_size;jj++) {u_val[jj]=0.0; } 
      for (k=0; k<K_size; k++) /* switch data for solve_eqs: Mu -> F */
      {
        mval = F_val[k]  ;
        F_val[k]  = Mu_val[k] ;
        Mu_val[k] = mval ;
      }
      for (k=0; k<K_len; k++) /* switch data for solve_eqs: M -> K */
      {
        mval     = K_val[k] ;
        K_val[k] = M_val[k] ;
        M_val[k] = mval ;
      }
      /* equation solver*/
      solve_eqs(K_val, u_val, F_val, K_size); /* equation solver*/
      for (k=0; k<K_size; k++) /* switch data for solve_eqs F -> Mu */
      {
        mval = F_val[k] ;
        F_val[k]  = Mu_val[k] ;
        Mu_val[k] = mval ;
      }
      for (k=0; k<K_len; k++) /* switch data for solve_eqs: K -> M */
      {
        mval     = K_val[k] ;
        K_val[k] = M_val[k] ;
        M_val[k] = mval ;
      }

      if (i > 0)
      {
        if ((fabs(omega - omega0)/omega) <= (pow(10,(-2.0*0.01))))
        {
          switch(sol_mode)
          {
            case 1:
              fprintf(stderr," Eigenvalue [%i]: %f (in iteration %i)\n",j,
              (fabs(omega)),i+1);
              break;
            case 2:
              fprintf(stderr," Eigenvalue [%i]: %f (in iteration %i)\n",j,
              sqrt(fabs(omega))/(2.0*3.141592),i+1);
              break;
          }
          rv = 0 ; /* we are converged */
          break ;
        }
      }
      omega0 = omega ;
    }
    if (rv != 0)
    {
      fprintf(stderr," Computation of eigenvalue [%i] failed\n",j);
      goto memFree;
    }
    /* prepare next Gram-Schmidt: */
    for (k=0; k<K_size; k++) { eig_val[j-1][k] = u_val[k] ; }
  }

  /* TODO: some savings of results... */

memFree:
  for (i=0; i<num_res; i++)
  {
    free(eig_val[i]); eig_val[i] = NULL ;
  }
  free(eig_val); eig_val = NULL ;
#endif
  return(rv) ;
}

#ifndef GR
int main(argc, argv)
int argc ;
char *argv[];
{
  FILE *fw = NULL ;
  FILE *fo = NULL ;
  FILE *fd = NULL ;
  FILE *fp = NULL ;
#ifdef UI
	char  name[16] ;
	int   i ;
#endif

  fprintf(stderr,"\nDDFOR 1.1.2: direct stiffness method solver for statics of 2D frames.\n");
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

  if (argc < 4)
  {
    /* no output */
    fd = NULL;
  }
  else
  {
    if ((argv[3][0] != '0')&&(argv[3][0] != '-'))
    {
      /* write to file */
      if ((fd=fopen(argv[3],"w")) == NULL)
      {
        fprintf(stderr,"Failed to open output file!\n");
        fd = NULL ;
      }
    } else {
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
    if ((argv[4][0] != '0')&&(argv[4][0] != '-'))
    {
      /* write to file */
      if ((fp=fopen(argv[4],"w")) == NULL)
      {
        fprintf(stderr,"Failed to open gfx file!\n");
        fp = NULL ;
      }
    } else {
      fp = NULL ;
    }
  }

  if (alloc_kf() != 0 )
  {
    free_data();
    return(-1);
  }

#ifdef LARGE
	if (sol_mode == 0) /* linear solver */
	{
#endif
  	stiff(0,0); 
  	disps_and_loads(0);

  	fprintf(stderr,"\nSolution: \n");
    solve_eqs(K_val, u_val, F_val, K_size); /* equation solver*/
  	fprintf(stderr,"End of solution. \n");

  	if (fo != NULL) 
  	{ 
    	results(fo);
#ifndef NO_PSEUDO_GFX
    	fprintf(fo,"\nScheme of structure with nodes numbers: \n\n");
    	pseudo_geom(fo, 0);
    	fprintf(fo,"\nElements numbers:\n\n");
    	pseudo_geom(fo, -1);
#endif
    	fclose(fo);
  	}

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
#ifdef LARGE
	}
	else /* eigenvalues solver */
	{
    if (sol_mode == 1) /* linear stability */
    {
  	  fprintf(stderr,"\nSolution (linear stability): \n");
  	  stiff(0,0); 
  	  disps_and_loads(0);
      solve_eqs(K_val, u_val, F_val, K_size); /* equation solver*/
      geom_stiff(0);
      for (i=0; i<K_len; i++) { K_val[i] = 0.0 ; }
      for (i=0; i<K_size; i++) 
      { 
        F_val[i] = 0.0 ; 
        u_val[i] = 1.0 ; 
      }
  	  stiff(0,0); 
  	  disps_and_loads(0);

      inv_iter(1) ; /* one is enough */
		  /* TODO */
  	  fprintf(stderr,"End of solution. \n");
    }
    else /* sol_mode == 2, modal analysis */
    {
  	  fprintf(stderr,"\nSolution (modal analysis): \n");
  	  stiff(0,0); 
      for (i=0; i<K_size; i++) { u_val[i] = 1.0 ; }
  	  disps_and_loads(0);
      inv_iter(2) ;
      /* TODO */
  	  fprintf(stderr,"End of solution. \n");
    }
	}
#endif

  free_sol_data();
  free_data();
  return(0);
}
#endif

/* end of ddfor.c */
