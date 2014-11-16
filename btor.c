/*
   File name: btor.c
   Date:      2014/11/16 20:08
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
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI    3.141592653589793238462643383279502884

/** Computes It for torsion, type: 0=circle, 1=rectangle */
double i_torsion(int type, double b, double h)
{
  switch (type)
  {
    case 0:
      return(PI/32.0*pow(b,4));
      break ;
    case 1:
    default:
      if ((h<=0)||(b<=0)) return(0.0);
      return((1.0/3.0)*(1.0-0.63*b/h*(1.0-pow(b,4)/(12*pow(h,4))))*pow(b,3)*h);
      break;
  }
  return(0.0);
}

int main (int argc, char *argv[])
{
  int type ;
  double b, h, it ;

  if (argc <=1)
  {
    fprintf(stderr,"Torsion momen (It) computation. Use %s type b h.\n",argv[0]);
    fprintf(stderr,"  Where: type is 0..circle, 1..rectangle \n");
    fprintf(stderr,"      -          ---    -                \n");
    fprintf(stderr,"     / \\        |   |   h               \n");
    fprintf(stderr,"     \\ /        |   |   |               \n");
    fprintf(stderr,"      -          ---    -                \n");
    fprintf(stderr,"    |-b-|       |-b-|                    \n");
    return(0);
  }

  if (argc > 1) { type = atoi(argv[1]) ; }
  if (argc > 2) { b = atof(argv[2]) ; }
  if (type > 0)
  {
    if (argc > 3) { h=atof(argv[3]) ; }
  }
  else h = 0.0 ;

  it = i_torsion(type, b, h);
  fprintf(stdout,"It = %e\n", it);
  
  return(0);
}

/* end of btor.c */
