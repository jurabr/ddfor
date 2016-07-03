/*
   File name: dgen.c
   Date:      2016/06/30 18:30
   Author:    Jiri Brozovsky

   Copyright (C) 2016 Jiri Brozovsky

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

   Model data generator for "dslab":
   

   +-----------+
   |           |
   |           |
   |     + F   | b, m
   |           |
   |           |
   +-----------+
        a, n

   Notes: 
   * one force (F) in the center
   * w=0 supports on all edges
   * n,m must be even numbers
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float F ;   /* force size */
float a, b ; /* slab dimensions */
int   n, m ; /* f.e. mesh divisions */

float E  = 20e9 ;
float nu = 0.2 ;
float t  = 0.1 ;
int   gr = 0 ;

/* reads data from given stream */
int read_data(FILE *fr)
{
  fscanf(fr,"%e %e %e %d %d", &F, &a, &b, &n, &m);

  if ((((int)(n/2.0))*2) != n)  n = (((int)(n/2.0))*2) ;
  if ((((int)(m/2.0))*2) != m)  m = (((int)(m/2.0))*2) ;

  if (a <= 0) {return(-1);}
  if (b <= 0) {return(-1);}
  if (fabs(F) == 1e-9) {return(-1);}
  return(0);
}

/* writes model to given stream */
int comp_out(FILE *fw)
{
  int i, j ;
  float x, y, dx, dy;
  int n1, n2, n3, n4 ;

  /* nodes: */
  fprintf(fw,"%i\n", (n+1)*(m+1));
  dx = a/(float)n ;
  dy = b/(float)m ;

  for (j=0; j<(m+1); j++)
  {
    y = (float)j*dy ;
    for (i=0; i<(n+1); i++)
    {
      x = (float)i*dx ;
      fprintf(fw,"%e %e\n",x,y); 
    }
  }

  /* elements: */
  fprintf(fw,"%i\n", (n)*(m));

  for (j=0; j<m; j++)
  {
    for (i=0; i<n; i++)
    {
      n1 = (j*(n+1))+i+1 ;
      n2 = n1 + 1   ;
      n4 = n1 +(n+1)  ;
      n3 = n4 + 1  ;
      fprintf(fw,"%i %i %i %i %e %e %e %i\n",n1,n2,n3,n4,E,nu,t,gr); 
    }
  }

	/* supports: */
  fprintf(fw,"%i\n", (n*2)+(m*2));
	/* up and down */
	for (i=0; i<n+1; i++) { fprintf(fw,"%i 1 0 %i\n",i+1,gr); }
	for (i=0; i<n+1; i++) { fprintf(fw,"%i 1 0 %i\n",(m*(n+1))+i+1,gr); }
	/* left: */
	for (i=1; i<m; i++) { fprintf(fw,"%i 1 0 %i\n",i*(n+1)+1,gr); }
	for (i=1; i<m; i++) { fprintf(fw,"%i 1 0 %i\n",(i+1)*(n+1),gr); }

	/* force in the middle : */
	n1 = (int)(((n+1)*(m+1))/2+1);
  fprintf(fw,"1\n%i 1 %e %i\n",n1,F,gr); 

  return(0);
}

int main(int argc, char *argv[])
{
  int rv = 0 ;

  if ((rv=read_data(stdin)) != 0) {return(rv);}
  return(comp_out(stdout));
}

/* end of dgen.c */
