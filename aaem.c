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

/** Testing routine for B3 */
void test_B3(void)
{
	int i ;

	for (i=1; i<=10; i++)
		fprintf(stdout,"%2i  %e\n",i*i,J_B3(i, (double)(i*i), 1.0/365.0, 30e9));
}

/** main routine: */
int main(int argc, char *argv[])
{
	test_B3();
	return(0);
}

/* end of aaem.c */
