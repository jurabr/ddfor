/*
   File name: bdesign.c
   Date:      2014/11/10 18:03
   Author:    Jiri Brozovsky

   Copyright (C) 2014 Jiri Brozovsky

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   in a file called COPYING along with this program; if not, write to
   the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
   02139, USA.

   Finds the smallest posible IPE cross-section for given internal forces.
   The material is expected to be a S235 steel (fy=0.9*235MPa).
   Possible IPE: 80..330

   NOTE: the code uses a lot of simplifications so it is 
         SUITABLE   FOR   PRELIMINARY    DESIGN    ONLY!
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int assess_IPE(double M, double V, double N)
{
  int i ;
  double h, A, I, W, S, Sh,f,e,a, fy;
  double sigma_top, sigma_bot, tau_top, tau_max ;
  
  fy = 235e6 * 0.9 ; /* design value of yield stress */
  
  for (i=0; i<12; i++)
  {
    h = 80 + i*20 ;
    if ((h > 240)&&(h < 300)) h = 270 ; /* fix fo real h sizes */
    if (h >= 300) h = 330 ; /* fix fo real h sizes */

    /* approximations: */
    A = h*h*0.000338738  + h*0.0836899 -1.60321 ;
    I = h*h*h*0.000475247 -0.0512346*h*h + 212.463 ;
    a = 0.0148498*h + 2.62096;  /* wall thickness */
    f = 0.0249428*h + 3.39449 ; /* thicknes of top part  */
    e = 0.465021*h + 8.23748  ; /* width of the top part */

    S = ((h/2.0-a/2.0)*f)/1e6 ; /* static moment S: UNUSED ATM */
    Sh = S + (pow(h/2-a,2)/2.0*a)/1e6 ; /* static moment for center (h/2) */
    
    h = h/1e3 ; /* mm  -> m  */
    A = A/1e4 ; /* cm2 -> m2 */ 
    I = I/1e8 ; /* cm4 -> m4 */
    W = I/(0.5*h) ;
    
#ifdef DEVEL
    fprintf(stdout," h=%e, A=%e, I=%e, W=%e\n", h, A, I, W);
#endif

    /* Stresses: */
    sigma_top = N/A - M/W ;
    sigma_bot = N/A + M/W ;
    tau_top   = V*S  / I*a ; /* UNUSED ATM */
    tau_max   = V*Sh / I*a ;

#ifdef DEVEL
    fprintf(stdout," fy=%e| Sxt=%e, Sxb=%e, Tmax=%e\n", fy*0.9,sigma_top, sigma_bot, tau_max);
#endif

    if ((fabs(sigma_top)<fy)&&(fabs(sigma_bot)<fy)&&(fabs(tau_max)<fy/sqrt(3)))
    { 
      fprintf(stdout,"Found IPE%i\n",(int)(h*1000));
      break;
    }
#ifdef DEVEL
    else
    {
      fprintf(stdout,"IPE%i...\n",(int)(h*1000));
    }
#endif
  }
  
  return(0);
}

int main(int argc, char *argv[])
{
  double M = 0 ;
  double V = 0 ;
  double N = 0 ;

  /* TODO: ddfor interface? */

  if (argc <= 1)
  {
    fprintf(stderr,"\nIPE beam designer. Use \"%s M V N\"! Units: [N,m].\n",
            argv[0]); exit(0);
  }
  if (argc > 1) M = atof(argv[1]);
  if (argc > 2) V = atof(argv[2]);
  if (argc > 3) N = atof(argv[3]);

  assess_IPE(M, V, N) ;

  return(0);
}

/* end of bdesign.c */
