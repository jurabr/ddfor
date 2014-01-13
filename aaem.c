/*
 * File name: aaem.c
 * Date:      2014/01/13 14:44
 * Author:    Jiri Brozovsky
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
double fi_AAEM(double epos, double t, double t_1, double E28, double E_t1)
{
  return(E_t1*J_B3(epos, t, t1, E28) - 1.0);
}

/** Computes creep function for AAEM */
double ksi_AAEM(double epos, double t, double t_1, double E28, double E_t1)
{
  return( (E_t1 / (E_t1-R_B3(epos,t,t1,E28))) -
          (1.0  / fi_AAEM(epos,t,t1,E_t1,E_t1)) );
}

/** Computes Age Adjusted Effective Modulus value */
double E_AAEM(double epos, double t, double t_1, double E28, double E_t1)
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

/** main routine: */
int main(int argc, char *argv[])
{
	test_B3();
	return(0);
}

/* end of aaem.c */
