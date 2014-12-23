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

   Finds the smallest posible number of steel rods for concrete rectangle.
   The material is expected to be:
     10425 steel (fy=420/1.15MPa)
     C16/20 concrete (fc=16/15 MPa, E=27GPa)

   NOTE: the code uses a lot of simplifications so it is 
         SUITABLE   FOR   PRELIMINARY    DESIGN    ONLY!
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int assess_bet(double M, double b, double h, double R)
{
  double c, z, A, d, n,  fy ;
  int nn; 
  
  fy = 420e6 / 1.15 ; /* design value of yield stress */
  c  = 10 ;
  if (R*1.5 > c) c = 1.5*R ;
  
  d = h - 0.001*(c+0.5*R) ;
  if (d < 0.002*R)
  {
    fprintf(stderr, "Height (h) too small, at least: %f\n",h+0.002*R);
    return(-1);
  }

  z = 0.85 * d ;
  A = 1e3*M / (z*fy) ; /* design area */
  n = A / (0.25*3.14*R*R*0.000001);

#if DEVEL
  fprintf(stderr,"A=%e, As1=%e (num=%e)\n",A, 0.25*3.14*R*R*0.000001, n);
#endif

  if ((double)((int)n) < n) { nn = (int)n+1 ; }
  else                      { nn = (int)n ; }
  n = (double)nn ;

  if (((3+2*(n-1)+n)*R)>(b*1000.0))
  {
    fprintf(stderr, "Width (b) too small, at least: %f\n",((3+2*(n-1)+n)*R)/1000.0);
    return(-1);
  }

  fprintf(stdout,"Number of diam. \"%i\" bars N = %i\n",(int)R,nn); /* result */
  
  return(0);
}

int main(int argc, char *argv[])
{
  double M = 0 ;
  double b = 0 ;
  double h = 0 ;
  double R = 0 ;

  /* TODO: ddfor interface? */

  if (argc <= 1)
  {
    fprintf(stderr,"\nRC beam designer. Use \"%s M[kNm] b[m] h[m] diam[mm]\"!\n",
            argv[0]); exit(0);
  }
  if (argc > 1) M = atof(argv[1]);
  if (argc > 2) b = atof(argv[2]);
  if (argc > 3) h = atof(argv[3]);
  if (argc > 4) R = atof(argv[4]); /* steel diameter */

  assess_bet(M, b, h, R) ;

  return(0);
}

/* end of bdesign.c */
