/*
 * File name: aaem.c
 * Date:      2014/01/13 14:44
 * Author:    Jiri Brozovsky
 *
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
 *
 * Description: AAEM  method + B3 model (Bazant et al) for DDFOR code
 *
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"


/* INPUT DATA B3: */
double E28 ;

/* INPUT DATA for AAEM: */
double t1, t ; /* studied times: in YEARS */

/* *************************************** */

/** External stuff for ddfor.c: */
extern int    sol_mode ; /* 0=statics, 1=stability, 2=modal, 3=aaem */
extern int    n_nodes ;
extern int    n_elems ;
extern int    n_disps ;
extern int    n_nfors ;
extern int    n_eload ;

extern float *x_i ;
extern float *y_i ;

extern float *E ;
extern float *A ;
extern float *I ;
extern float *rho ;
extern int   *n1 ;
extern int   *n2 ;
extern int   *type ;

extern int   *d_n ; /* node */
extern int   *d_d ; /* direction 1=x 2=y 3=rot */
extern float *d_v ; /* size */

extern int   *f_n ; /* node */
extern int   *f_d ; /* direction 1=fx 2=fy 3=m */
extern float *f_v ; /* size */

extern int   *l_e ; /* node */
extern int   *l_d ; /* direction 1=x 2=y, 3=x global, 4=y global */
extern float *l_v1 ; /* size at beginning */
extern float *l_v2 ; /* size at end */

/* solution variables: */
extern double  *ke[] ;
extern double  *keg[] ;
extern double  *T[] ;
extern double  fe[];
extern double  feg[];
extern double  ueg[];
extern double  ue[];

extern int    K_len    ;
extern int   *K_sizes  ; 
extern int   *K_from   ;
extern int   *K_cols   ;
extern double *K_val   ;
extern double *F_val   ;
extern double *u_val   ;
extern double *Fr_val  ;
extern double *uu_val  ;
extern double *EE      ; 
extern double *E1      ; 

/** Reading of data from file */
extern int read_data(FILE *fw);
extern void free_sol_data();
extern int alloc_kf();
extern void free_data();
extern void stiff(int eg, int lc);
extern void disps_and_loads(int lc);
extern int solve_eqs(); /* conjugate gradient solver */
extern void results(FILE *fw);
extern void eint_results(FILE *fw);
extern void gfx_results(FILE *fw);


/** Computes compliance function for B3 (no drying) */
double J_B3(int epos/*unused*/, double t, double t1, double E28)
{
	double n,m, alpha, psi ; /* inputs        */
	double E0, qs, Cd;       /* computed data */

	/* predefined default data (Bazant, Jirasek): */
	n     = 0.1 ;
	m     = 0.5 ;
	alpha = 0.001 ;
	psi   = 0.3 ;

	/* computation of other parameters: */
	E0 = E28 / 0.6  ;
	qs = 11.4 / E28 ;
	Cd = 0.0 ; /* drying is not used at the moment */

	/* J(t,t1) value: */
	return(1.0/E0+qs*log(1.0+psi*((pow(t1,-m)+alpha)*pow(t-t1,n)))+Cd);
}

/** Computes approximation of R(t,t1) for concrete */
double R_B3(int epos/*unused*/, double t, double t1, double E28)
{
  double dt ;

  dt = t - t1 ; /* time difference */

  return(0.992/J_B3(epos,t,t1,E28)-((0.115/J_B3(epos,t,t-1.0,E28))  
         *((J_B3(epos,dt,t1,E28)/J_B3(epos,t,dt,E28))-1.0)) );
}

/** Computes aging function for AAEM */
double fi_AAEM(int epos, double t, double t_1, double E28, double E_t1)
{
  return(E_t1*J_B3(epos, t, t1, E28) - 1.0);
}

/** Computes creep function for AAEM */
double ksi_AAEM(int epos, double t, double t_1, double E28, double E_t1)
{
  return( (E_t1 / (E_t1-R_B3(epos,t,t1,E28))) -
          (1.0  / fi_AAEM(epos,t,t1,E_t1,E_t1)) );
}

/** Computes Age Adjusted Effective Modulus value */
double E_AAEM(int epos, double t, double t_1, double E28, double E_t1)
{
  return(
    E_t1/( 1.0+ksi_AAEM(epos,t,t1,E_t1,E_t1)*fi_AAEM(epos,t,t1,E_t1,E_t1)));
}

/** TESTING ROUTINES ----------------------------- */

/** Testing routine for B3 */
void test_B3(void)
{
	int i  ;
  double days ;

	for (i=1; i<=10; i++)
  {
    days = (double)(i*365) ; /* 1 year */
		fprintf(stdout,"%3.0f %3.10e %3.10e\n",days,
        J_B3(i, days, 28.0, 30e9), R_B3(i, days, 28.0, 30e9));
  }
}

/** ANALYSIS OF FRAMES --------------------------- */

/** AAEM for 2D frames (uses DDFOR) */
int aaem_frame(int argc, char *argv[])
{
  FILE *fw = NULL ;
  FILE *fo = NULL ;
  FILE *fd = NULL ;
  FILE *fp = NULL ;
  int   i ;
  double t, t1 ;

  /* TODO: replace testing data with user-changeable stuff: */
  t1 = 28.0 ;      /* 28 days */
  t  = 5.0*365.0 ; /* 5 years */

  sol_mode = 3 ; /* AAEM: needed for data allocations */

  fprintf(stderr,"\nDDFOR/AAEM 0.1: time-dependent analysis of 2D concrete frames.\n");
  fprintf(stderr,"  See for details: http://github.com/jurabr/ddfor\n\n");

  if (argc < 2)
  {
    /* standard input */
    fprintf(stderr,"\nInteractive input not possible!\n");
    fprintf(stderr,"Use: %s data_file!\n\n",argv[0]);
    return(-1);
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

  /* SOLUTION START -------------- */
  /* solution for initial time t1: */
  fprintf(stderr,"\nSolution: \n");
  for (i=0; i<n_elems; i++) 
  { 
    EE[i] = E[i] ; /* save initial E28 */
    E[i] = E_AAEM(i, t, t1, EE[i], EE[i]) ; /* E1 */
    E1[i] = E[i] ; /* save E1 */
  } /* TODO: second EE[i] replace with time-dependent  function! */
  
  /* first AAEM run: */
 	stiff(0,0); 
  disps_and_loads(0); /* only initial loads are here */
  solve_eqs(); 
  for (i=0; i<K_len; i++) { uu_val[i] = u_val[i] ; } /* results(t1) */

  for (i=0; i<K_len; i++) /* clean previous data*/
  {
    K_val[i]  = 0.0 ;
    F_val[i]  = 0.0 ;
    u_val[i]  = 0.0 ;
  }
  
  /* second AAEM run (for time "t") - 1st part */
  for (i=0; i<n_elems; i++) /* prepare E" for given time */
  {
    E[i] = E_AAEM(i, t, t1, EE[i], E1[i]) ; /* E(t) */
  } 
  stiff(0,1); 
  disps_and_loads(1); /* TODO: use t1 loads here!!! */
  solve_eqs(); 
  for (i=0; i<K_len; i++) { Fr_val[i] = u_val[i] ; } /* partial results(t) */

  /* second AAEM run (for time "t") - 2nd part */
  for (i=0; i<n_elems; i++) /* prepare E" */
  {
    E[i] = E1[i] / fi_AAEM(i,t,t1,EE[i],E1[i]) ;
  } 
  stiff(0,1); 
  disps_and_loads(1);  /* TODO: use t1 loads here!!! */
  solve_eqs(); 
  for (i=0; i<K_len; i++) { u_val[i] += Fr_val[i] ; } /* final results */

  fprintf(stderr,"End of solution. \n");
  /* SOLUTION END ---------------- */

    /* Outputs: */
  	if (fo != NULL) 
  	{ 
    	results(fo);
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

  free_sol_data();
  free_data();
  return(0);
}

/** main routine: */
int main(int argc, char *argv[])
{
	test_B3();
	return(0);
}

/* end of aaem.c */
