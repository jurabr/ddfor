#include <stdio.h>
#include <stdlib.h>
#include <graphics.h>
#include <bios.h>

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
double xcorr= 1.0 ; /* x size correction */
int gr_size = 3 ; /* size of symbols */

double g_dx = 0.0 ;
double g_dy = 0.0 ;
double gmul = 1.0 ;

/** Reading of data from file */
extern int read_data(FILE *fw);

/** free K,F */
extern void free_sol_data();

/* Allocates space for linear system */
extern int alloc_kf();

extern double norm_K();

extern double vec_norm(double *a, int len);

extern int solve_eqs();

/** Local stiffness matrix */
extern void stiff_loc(int type, float E, float A, float I, float l);

/* set transformation matrix to zero */
extern void tran_zero();

/** Transformation matrix */
extern void tran(double s, double c);

/* fils "ke" with content of "keg" */
extern void ke_switch();

extern void ke_to_keg(double s, double c);

/* puts data to the right place in K */
extern void md_K_add(int row, int col, float val);

/* loads on nodes */
extern void one_eload(int epos, double na, double nb, double va, double vb, double L);

/* computes stiffness matrix of the structure */
extern void stiff();

/* supports */
extern void add_one_disp(int node, int dir, float val);

/* forces */
extern void add_one_force(int node, int dir, float val);

extern void disps_and_loads();

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
extern void in_gfx(int type, int epos, int div, int ppos, double *mult, double *vx, double *vy);


/* sets GFX mode */
void gr_init(void)
{
  setvmode(6); /* CGA 640x200 BW */
  okraj =   15 ;
  imaxx  = 320 ; /* limited by a reason - real screen parameters */
  imaxy  = 200 ;
  xcorr  = 2.0 ; /* for MC600, 200LX */
  gr_size= 3   ;
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
  
  if (muly == 0.0) muly = mulx ;
  if (mulx == 0.0) mulx = muly ;

  if (mulx < muly) { gmul = mulx ; }
  else             { gmul = muly ; }

  /*fprintf(stderr,"GMUL (%f %f): %f /ix=%f iy=%f\n",mulx,muly,gmul,ix,iy); */
}

/** Computation of screen coordinates */
int x_pos(double x) { return((int)(xcorr*(x - g_dx)*gmul))+okraj; }
int y_pos(double y) { return(imaxy - ((int)((y - g_dy)*gmul))-okraj); }

/** plots supports */
void plot_disp(int node, int type, double val, int size)
{
  int x,y;

  x = x_pos(x_i[node-1]);
  y = y_pos(y_i[node-1]);

  switch (type)
  {
    case 1: move_to(x,y); line_to(x-2*size*xcorr, y);
            move_to(x-2*size*xcorr, y+size);
            line_to(x-2*size*xcorr, y-size);
            break ;
    case 2: move_to(x,y); line_to(x,y+2*size);
            move_to(x+size*xcorr, y+2*size);
            line_to(x-size*xcorr, y+2*size);
            break ;
    case 3: move_to(x-size*xcorr,y-size);
            line_to(x+size*xcorr,y-size);
            line_to(x+size*xcorr,y+size);
            line_to(x-size*xcorr,y+size);
            line_to(x-size*xcorr,y-size);
            break ;
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
}

/** prints internal forces in elements gfx-friendly form */
void gfx_plot_results(FILE *fw)
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


/** main routine */
int main(int argc, char *argv[])
{
  FILE *fw = NULL ;
  FILE *fo = NULL ;
  FILE *fd = NULL ;

  fprintf(stderr,"\nDDFOR 1.0.4: direct stiffness method solver for statics of 2D frames.\n");
  fprintf(stderr,"  See for details: http://github.com/jurabr/ddfor\n\n");

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

  if (argc < 3)
  {
    /* standard output */
    fo = stdout;
  }
  else
  {
    /* write to file */
    if ((fo=fopen(argv[2],"w")) == NULL)
    {
      fprintf(stderr,"Failed to open output file!\n");
      fo = stdout ;
    }
  }

  if (argc < 4)
  {
    /* no output */
    fd = NULL;
  }
  else
  {
    /* write to file */
    if ((fd=fopen(argv[3],"w")) == NULL)
    {
      fprintf(stderr,"Failed to open output file!\n");
      fd = NULL ;
    }
  }

  if  (alloc_kf() != 0 )
  {
    free_data();
    return(-1);
  }

  stiff(); 
  disps_and_loads();

  fprintf(stderr,"\nSolution: \n");
  solve_eqs();
  fprintf(stderr,"End of solution. \n");

  /* results(fo); */
  fclose(fo);

  if (fd != NULL) 
  { 
    /* eint_results(fd);  */
    fclose(fd);
  }

  /** Graphics: */
  gr_init();
  set_minmax();
  plot_struct();
  getchar();
  gr_stop();

  free_sol_data();
  free_data();
  return(0);
}
