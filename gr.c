#include <stdio.h>
#include <stdlib.h>
#include <graphics.h>
#include <bios.h>
#include <math.h>

/** External stuff for ddfor.c: */
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
extern int   *K_sizes  ; /* lenghts of K's rows */
extern int   *K_from   ;
extern int   *K_cols   ;
extern double *K_val   ;
extern double *F_val   ;
extern double *u_val   ;

extern double *M    ;
extern double *r    ;
extern double *z    ;
extern double *p    ;
extern double *q    ;


/* GFX parameters: */
int okraj   = 0 ; /* border */
int imaxx   = 0 ;
int imaxy   = 0 ;
int gr_size = 3 ; /* size of symbols */

double g_dx = 0.0 ;
double g_dy = 0.0 ;
double gmul = 1.0 ;
double gsiz = 0.0 ;
double grat = 1.0 ; /* g_x/g_y ratio */
double zoom = 1.0 ; /* zoom          */
double fmax = 1.0 ; /* max. size of loads */
double efmax = 1.0 ; /* max. size of loads */

/** Reading of data from file */
extern int read_data(FILE *fw);

/** free K,F */
extern void free_sol_data();

/* Allocates space for linear system */
extern int alloc_kf();

extern double norm_K(); /* matrix norm */

extern double vec_norm(double *a, int len); /* vector norm */

extern int solve_eqs(); /* conjugate gradient solver */

/** Local stiffness matrix */
extern void stiff_loc(int type, float E, float A, float I, float l);

/* set transformation matrix to zero */
extern void tran_zero();

/** Transformation matrix */
extern void tran(double s, double c);

/* fills "ke" with content of "keg" */
extern void ke_switch();

extern void ke_to_keg(double s, double c); /* local->global element matrix */

/* puts data to the right place in K */
extern void md_K_add(int row, int col, float val);

/* loads on nodes */
extern void one_eload(int epos, double na, double nb, double va, double vb, double L);

/* computes stiffness matrix of the structure */
extern void stiff(int eg, int lc);

/* supports */
extern void add_one_disp(int node, int dir, float val);

/* forces */
extern void add_one_force(int node, int dir, float val);

extern void disps_and_loads(int lc);

/** Frees all allocated data */
extern void free_data();

extern void u_to_ue(double s, double c);

/** Find all element load in give direction */
extern void get_eloads(int epos, int dir, double *na, double *nb);

/** Compute internal force (N, V, M) for given point of beam */
extern double in_force(int type, int epos, int div, int ppos);

/** Local results */
extern void res_loc(int epos);

/** computes and prints results */
extern void results(FILE *fw);

/** computes and prints internal forces in elements */
extern void eint_results(FILE *fw);

/** Find (approximation of) extreme values of N, V, M on beam */
extern int beam_max(int type, int epos, int div, double *nmax, double *npos, double *vmax, double *vpos, double *mmax, double *mpos);

/** Compute internal force (N, V, M) for given point of beam
 *  This routine is meant for plotting tools */
extern void in_gfx(int type, int epos, int div, int ppos, double mult, double *vx, double *vy);

void plint(int num) /* plot integer value */
{
  char s[10] ;
  int  i ;

  if (num > 9e9) return;
  for (i=0; i<10; i++) s[i]='\0';
  sprintf(s,"%i",num);
  plots(s);
}

void pldbl(double num) /* plot double value */
{
  char s[10] ;
  int  i ;

  if (abs(num) > 9e7) return;
  for (i=0; i<10; i++) s[i]='\0';
  sprintf(s,"%.0f",num);
  plots(s);
}

/* sets GFX mode */
void gr_init(void)
{
  setvmode(5); /* CGA 320x200 BW */
  okraj =   15 ;
  imaxx  = 320 ; 
  imaxy  = 200 ;
  gr_size= 3   ;
  grat = (double)imaxx / (double)imaxy ;
}

/** restore text mode */
void gr_stop(void)
{
  setvmode(DEFAULTMODE); /* restore original video mode */
  return(0);
}

/* gets dimension of the structure and the coefficients */
void set_minmax(void)
{
  double minx, maxx, miny, maxy ;
  double lx, ly, mulx, muly, ix, iy ;
  int i ;
  int cent = 0 ;
  double cmul = 0.0 ;

  maxx = x_i[0]; maxy = y_i[0];
  minx = x_i[0]; miny = y_i[0];

  for (i=0; i<n_nodes; i++) /* get limits */
  {
    if (x_i[i] > maxx) { maxx = x_i[i] ; }
    if (x_i[i] < minx) { minx = x_i[i] ; }
    if (y_i[i] > maxy) { maxy = y_i[i] ; }
    if (y_i[i] < miny) { miny = y_i[i] ; }
  }
  g_dx = (-1.0)*minx ; /* starting points */
  g_dy = (-1.0)*miny ;

  lx = abs(maxx - minx) ; /* multipliers */
  ly = abs(maxy - miny) ;
  mulx = 0.0 ; muly = 0.0 ;
  ix = (double)(imaxx - (2*okraj)); 
  iy = (double)(imaxy - (3*okraj)); 
  if (lx > 0.0) { mulx = ix/lx; }
  if (ly > 0.0) { muly = iy/ly; }
  
  if (muly<=0.3*mulx) 
  {
    cent=2;
    if (muly!=0.0) cmul=0.8*mulx/muly;
    else cmul=0.3;
  }
  if (mulx<=0.3*muly)
  {
    cent=1;
    if (mulx!=0.0) cmul=0.8*muly/mulx;
    else cmul=0.3;
  }

  if (muly == 0.0) muly = mulx ;
  if (mulx == 0.0) mulx = muly ;

  if (mulx < muly) { gmul = mulx ; } /* gfx multiplier */
  else             { gmul = muly ; }

  if (lx > ly) { gsiz = lx ; } /* size of structure */ 
  else         { gsiz = ly ; }
  
  /* adjust position if lx or ly is very small */
  switch(cent)
  {
    case 1: g_dx -= (cmul*gsiz) ; break ;
    case 2: g_dy -= (cmul*gsiz) ; break ;
  }
}

/** Computation of screen coordinates */
int x_pos(double x) { return((int)((x + g_dx)*gmul*zoom))+okraj; }
int y_pos(double y) { return(imaxy - ((int)((y - g_dy)*gmul*zoom))-okraj); }

/** plots supports */
void plot_disp(int node, int type, double val, int size)
{
  int x,y;

  x = x_pos(x_i[node-1]);
  y = y_pos(y_i[node-1]);

  switch (type)
  {
    case 1: move_to(x,y); line_to(x-2*size, y);
            move_to(x-2*size, y+size);
            line_to(x-2*size, y-size);
            break ;
    case 2: move_to(x,y); line_to(x,y+2*size);
            move_to(x+size, y+2*size);
            line_to(x-size, y+2*size);
            break ;
    case 3: move_to(x-size,y-size);
            line_to(x+size,y-size);
            line_to(x+size,y+size);
            line_to(x-size,y+size);
            line_to(x-size,y-size);
            break ;
  }
}

/** plots one force  */
void plot_force(int node, int type, double val, int size0)
{
  int x,y;
  int size = 10 ;

  size = (int)(abs((val)/fmax) * (double)size0) ;
  if (size < 6) {size = 6;}

  x = x_pos(x_i[node-1]);
  y = y_pos(y_i[node-1]);

  switch (type)
  {
    case 1: /* horizontal force */
            if (val > 0.0)
            {
              move_to(x+4*size, y);
              line_to(x,y);
              line_to(x+size, y-size);
              move_to(x+size, y+size);
              line_to(x, y);

              move_to(x+2*size, y-2*size);
              pldbl(val);
            }
            else
            {
              move_to(x-4*size, y);
              line_to(x,y);
              line_to(x-size, y-size);
              move_to(x-size, y+size);
              line_to(x, y);

              move_to(x-4*size, y-2*size);
              pldbl(val);
            }
            break ;
    case 2: /* vertical force */
            if (val > 0.0)
            {
              move_to(x, y+4*size);
              line_to(x,y);
              line_to(x-size, y+size);
              move_to(x+size, y+size);
              line_to(x, y);

              move_to(x, y+4*size);
              pldbl(val);
            }
            else
            {
              move_to(x, y-4*size);
              line_to(x,y);
              line_to(x-size, y-size);
              move_to(x+size, y-size);
              line_to(x, y);

              move_to(x, y-5*size);
              pldbl(val);
            }
            break ;
    case 3: /* bending moment */
            if (val > 0.0)
            {
              move_to(x-2*size, y);
              line_to(x-2*size, y-2*size);
              line_to(x-1*size, y-3*size);
              line_to(x+1*size, y-3*size);
              line_to(x+2*size, y-2*size);
              line_to(x+2*size, y);
              move_to(x+1*size, y-size);
              line_to(x+2*size, y);
              line_to(x+3*size, y-size);

              move_to(x-size, y-4*size);
              pldbl(val);
            }
            else
            {
              move_to(x+2*size, y);
              line_to(x+2*size, y-2*size);
              line_to(x+1*size, y-3*size);
              line_to(x-1*size, y-3*size);
              line_to(x-2*size, y-2*size);
              line_to(x-2*size, y);
              move_to(x-1*size, y-size);
              line_to(x-2*size, y);
              line_to(x-3*size, y-size);

              move_to(x-size, y-5*size);
              pldbl(val);
            }
            break ;
  }
}
/** plot one element load */
void plot_eload (int elem, int type, double v1, double v2, int size)
{
  double x1,x2,y1,y2,c,s ;
  double L, lenx, vx, vy, vx0, vy0 ;
  double val = 0.0 ;
  double mult = 5.0 ;
  int div = 10 ;
  int i ;

  x1 = (double)x_pos(x_i[n1[elem]-1]) ;
  y1 = (double)y_pos(y_i[n1[elem]-1]) ;
  x2 = (double)x_pos(x_i[n2[elem]-1]) ;
  y2 = (double)y_pos(y_i[n2[elem]-1]) ;

  L = sqrt ( ((y2-y1)*(y2-y1)) + ((x2-x1)*(x2-x1)) ) ;
  s = (y2-y1)/L ;
  c = (x2-x1)/L ;

  switch(type)
  {
    case 1: /* TODO  local x */ 
      break ;
    case 2: /* TODO  local y */
        move_to(x1,y1);
        for (i=0; i<=div; i++)
        {
          lenx  = L*((double)((double)i/(double)(div))) ;

          vx0 = x1 + lenx*c ;
          vy0 = y1 + lenx*s ;

          val =  v1+(v2-v1)*((double)((double)i/(double)(div))) ;
          mult = (val/efmax)*(double)size ;
          vx = vx0 - mult*s ;
          vy = vy0 - mult*c ;
          line_to((int)vx,(int)vy);
          line_to((int)vx0,(int)vy0);
          move_to((int)vx,(int)vy);
        }
      break ;
    case 3: /* TODO  global x */ break ;
    case 4: /* TODO  global y */ break ;
    default: break ;
  }
}

/** plot nodes */
void plot_node(double x_i, double y_i)
{
  int x,y;

  x = x_pos(x_i);
  y = y_pos(y_i);

  move_to(x,y-2);
  line_to(x+2,y);
  line_to(x,y+2);
  line_to(x-2,y);
  line_to(x,y-2);
}

/** plot loads */
void plot_forces(void)
{
  int i ;

  /* multiplier for forces: */
  fmax = 0 ;
  efmax = 0 ;
  for (i=0; i<n_nfors; i++)
    { if (abs(f_v[i]) > fmax) {fmax = abs(f_v[i]);} }
  for (i=0; i<n_eload; i++)
  { 
    if (abs(l_v1[i]) > efmax) {efmax = abs(l_v1[i]);} 
    if (abs(l_v2[i]) > efmax) {efmax = abs(l_v2[i]);} 
  }

  for (i=0; i<n_nfors; i++)
  {
    if (f_v[i] != 0.0) plot_force(f_n[i], f_d[i], f_v[i], (int)(0.5*okraj)) ;
  }
  for (i=0; i<n_eload; i++)
  {
    if ((l_v1[i] != 0.0) || (l_v2[i] != 0.0))
      plot_eload(l_e[i]-1, l_d[i], l_v1[i], l_v2[i], (int)(0.5*okraj)) ;
  }
}

/** plot geometry */
void plot_struct(void)
{
  int i,j ;

  /* nodes: */
  for (i=0; i<n_nodes; i++) { plot_node(x_i[i],y_i[i]); }
  
  /* elements: */
  for (i=0; i<n_elems; i++)
  {
    move_to(x_pos(x_i[n1[i]-1]), y_pos(y_i[n1[i]-1]));
    line_to(x_pos(x_i[n2[i]-1]), y_pos(y_i[n2[i]-1]));
  }
  
  /* supports: */
  for (i=0; i<n_disps; i++) { plot_disp(d_n[i],d_d[i],d_v[i],gr_size); }

  /* forces */
  plot_forces();
}

/** plot element (with numbers) */
void gfx_plot_elements(void)
{
  int i ;

  for (i=0; i<n_elems; i++)
  {
    move_to(x_pos(x_i[n1[i]-1]), y_pos(y_i[n1[i]-1]));
    line_to(x_pos(x_i[n2[i]-1]), y_pos(y_i[n2[i]-1]));

    move_to(x_pos(0.5*(x_i[n1[i]-1]+x_i[n2[i]-1])), 
            y_pos(0.5*(y_i[n1[i]-1]+y_i[n2[i]-1]))-gr_size);
    plint(i+1);
  }
}

/** plot nodes (with numbers) */
void gfx_plot_nodes(void)
{
  int i,j ;

  /* nodes: */
  for (i=0; i<n_nodes; i++) 
  { 
    plot_node(x_i[i],y_i[i]); 
    move_to(x_pos(x_i[i])+gr_size,y_pos(y_i[i])-gr_size);
    plint(i+1);
  }
}

/** prints internal forces in elements gfx-friendly form */
void gfx_plot_results(int type)
{
  int i,j ;
  int div = 10 ;
  double x,y, mx,my ;
  double mult = 0.0;
  double fmult = 0.0;
#if 0
  double dmult = 0.0;

  /* multiplier for deformations */
  dmult = abs(u_val[0]) ;
  for (i=1; i<K_len; i++)
    { if (dmult < abs(u_val[i])) { dmult = abs(u_val[i]); } }
  if (dmult > 1e-6) { dmult = 0.3 * (gsiz/dmult) ;  } 
#endif

  mx = 0.0 ;
  for (i=0; i<n_elems; i++) /* compute result multiplier */
  {
    res_loc(i);
    for (j=0; j<=div; j++)
    {
      mx = abs(in_force(type, i, div, j));
      if (mx > fmult) { fmult = mx ; }
    }
  }
  if (fmult != 0.0) 
  { 
    mult = ((double)okraj/(double)imaxx)*(1.0/grat)*0.9*(gsiz/fmult); 
  }
  else { mult = 0.0 ; }

  for (i=0; i<n_elems; i++) /* plot results */
  {
    res_loc(i);
    for (j=0; j<=div; j++)
    {
      in_gfx(0, i, div, j, mult, &x,  &y);
      in_gfx(type, i, div, j, mult, &mx,  &my);
      if (j == 0 ) 
      { 
        move_to(x_pos(x), y_pos(y)); 
        line_to(x_pos(mx), y_pos(my)); 
      }
      else
      {
        line_to(x_pos(mx), y_pos(my)); 
        if (type != 4) {
        line_to(x_pos(x), y_pos(y));}   /* these 2 lines are for filling */
        move_to(x_pos(mx), y_pos(my)); 
      }
    }
  }
}

/** main routine */
int main(int argc, char *argv[])
{
  FILE *fw = NULL ;
  int   c  = 0 ;

  fprintf(stderr,"\nDDFOR/GFX 0.2.4: direct stiffness method solver for statics of 2D frames.\n");
  fprintf(stderr,"  See for details: http://github.com/jurabr/ddfor\n\n");

  /** normal program run: */
  if (argc < 2)
  {
    fprintf(stderr,"No input data file!\n");
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
      read_data(fw);
      fclose(fw);
    }
  }

  if  (alloc_kf() != 0 )
  {
    free_data();
    return(-1);
  }

  stiff(0,0); 
  disps_and_loads(0);

  fprintf(stderr,"\nSolution: \n");
  solve_eqs();
  fprintf(stderr,"End of solution. \n");

  /** Graphics: */
  gr_init();
  set_minmax();
  plot_struct();
  /* keyboard control loop for gfx: */
  while (c != 27)
  {
    c = getch();
    switch(c)
    {
      case 102: /* F,f: forces, loads */
      case  70: /* TODO */
        break;
      case 43: /* zoom += */
      case 61: 
        zoom += 0.1 ; if (zoom>10) zoom = 10 ;
        clrscrn2(0);
        plot_struct(); 
        break ;
      case 45: /* zoom -  */
        zoom -= 0.1 ; if (zoom<=0.0) zoom = 0.1 ;
        clrscrn2(0);
        plot_struct(); 
        break ;
      case 110: 
      case  78:
        gfx_plot_nodes(); 
        break;
      case  97: /* normal forces */  
      case  65:
        clrscrn2(0);
        gfx_plot_results(0); 
        gfx_plot_results(1); 
        break ;
      case 118: /* shear forces */ 
      case  86:
        clrscrn2(0);
        gfx_plot_results(0); 
        gfx_plot_results(2); 
        break ; 
      case  77: /* bending moments */ 
      case 109:
        clrscrn2(0);
        gfx_plot_results(0); 
        gfx_plot_results(3);
        break ;
      case 119: /* deformations */
      case  87:
        clrscrn2(0);
        /*gfx_plot_results(0); */
        gfx_plot_results(4);
        break ;
      case 101:/* elements */
      case  69:
        gfx_plot_elements();
        break ;
      case 113: /* quit */
      case  27:
        c = 27 ; break ;
      case 105: /* structure */
      case  73:
        clrscrn2(0);
        plot_struct(); 
        break ; 
      case  83:
      case 115:
      default: 
        clrscrn2(0);
        gfx_plot_results(0); 
        break ; 
    }
  }
  gr_stop();

  free_sol_data();
  free_data();
  return(0);
}
