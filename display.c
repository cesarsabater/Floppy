/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <math.h>

/* macros */

#define IX(i,j) ((i)+(N+2)*(j))
#define max(x,y) (((x) > (y)) ? (x) : (y))

/* external definitions (from solver.c) */

typedef void (*stepfun)(int,float**,float**,float**,float**,float,float);

extern stepfun get_dens_step();
extern stepfun get_vel_step();
extern stepfun get_dens_step_orig();
extern stepfun get_vel_step_orig();
extern void alloc_data();
extern void refresh_grid(float**);
extern void test();
extern int **grid;
extern int cflag;
extern int G;
extern int slot_size;

/* global variables */

stepfun vel_step; 
stepfun dens_step;
stepfun dens_step_orig;
stepfun vel_step_orig;

int N;
static float dt, diff, visc;
static float force, source;
static int dvel;
int pause;
int iter;

// optimized grid
float ** u, ** v, ** u_prev, ** v_prev;
float ** dens, ** dens_prev;
float ** vel_max, ** dens_max;

// original grid
float ** uo, ** vo, ** uo_prev, ** vo_prev;
float ** denso, ** denso_prev;

static int win_id;
static int win_x, win_y;
static int omx, omy, mx, my;


/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

static void free_matrix(float **m) {
	int i;
	for (i = 0; i < N+2; i++) 
		if (m[i]) 
			free(m[i]);
	free(m);
}

static void free_data (void)
{
	if ( u ) free_matrix ( u );
	if ( v ) free_matrix ( v );
	if ( u_prev ) free_matrix ( u_prev );
	if ( v_prev ) free_matrix ( v_prev );
	if ( dens ) free_matrix ( dens );
	if ( dens_prev ) free_matrix ( dens_prev );
	if (dens_max) free_matrix (dens_max); 
	if (vel_max) free_matrix(vel_max);
	if ( uo) free_matrix ( uo );
	if ( vo ) free_matrix ( vo );
	if ( uo_prev ) free_matrix ( uo_prev );
	if ( vo_prev ) free_matrix ( vo_prev );
	if ( denso ) free_matrix ( denso );
	if ( denso_prev ) free_matrix ( denso_prev );
}

static void clear_data (void)
{
	int i, j, size=N+2;
	for ( i=0 ; i<size ; i++ ) for (j=0 ; j<size; j++) {
		dens_max[i][j] = vel_max[i][j] = u[i][j] = v[i][j] = u_prev[i][j] = v_prev[i][j] = dens[i][j] = dens_prev[i][j] = 0.0f;
		uo[i][j] = vo[i][j] = uo_prev[i][j] = vo_prev[i][j] = denso[i][j] = denso_prev[i][j] = 0.0f;
	}
}

static float ** alloc_matrix() { 
	int i, size=N+2;
	float **m; 
	
	m = (float**) malloc ( size * sizeof(float*)); 

	if ( !m ) {
		fprintf ( stderr, "cannot allocate data m\n" );
		return NULL;
	}

	for (i = 0; i < size; i++)  { 
		m[i] = (float*) malloc(size * sizeof(float));
		if (!(m[i])) {
			fprintf ( stderr, "cannot allocate data\n" );
			return NULL; 
		}
	}

	return m; 
}

static int allocate_data ( void )
{
	// optim
	u				  = alloc_matrix(); v 			  = alloc_matrix();
	u_prev	  = alloc_matrix(); v_prev	  = alloc_matrix();
	dens		  = alloc_matrix(); dens_prev	= alloc_matrix();
	// some statistic values
	vel_max   = alloc_matrix();
	dens_max  = alloc_matrix();
	// orig
	uo				  = alloc_matrix(); vo 			  = alloc_matrix();
	uo_prev	  = alloc_matrix(); vo_prev	  = alloc_matrix();
	denso		  = alloc_matrix(); denso_prev	= alloc_matrix();
	
	if ( !u || !v || !u_prev || !v_prev ||!vel_max || !dens || !dens_prev || !dens_max
		|| !uo || !vo || !uo_prev || !vo_prev || !denso || !denso_prev) {
		fprintf ( stderr, "cannot allocate data\n" );
		return 0;
	}
	
	return 1;
}


/*
  ----------------------------------------------------------------------
   OpenGL specific drawing routines
  ----------------------------------------------------------------------
*/

static void pre_display ( void )
{
	glViewport ( 0, 0, win_x, win_y );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
	glutSwapBuffers ();
}

static void draw_velocity ( void )
{
	int i, j;
	float x, y, h;

	h = 1.0f/N;

	glColor3f ( 1.0f, 1.0f, 1.0f );
	glLineWidth ( 1.0f );

	glBegin ( GL_LINES );

		for ( i=1 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=1 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				glVertex2f ( x, y );
				glVertex2f ( x+u[i][j], y+v[i][j] );
			}
		}

	glEnd ();
}

static void draw_map(float **v)
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;

	h = 1.0f/N;

	glBegin ( GL_QUADS );

		for ( i=0 ; i<=N ; i++ ) {
			x = (i-0.5f)*h;
			for ( j=0 ; j<=N ; j++ ) {
				y = (j-0.5f)*h;

				d00 = v[i][j];
				d01 = v[i][j+1];
				d10 = v[i+1][j];
				d11 = v[i+1][j+1];

				glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
				glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
				glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
				glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
			}
		}
	
	glEnd ();
}

static void draw_density_orig ( void )
{
	draw_map(denso);
}

static void draw_density ( void )
{
	draw_map(dens);
}
/*
static void draw_max_dens(void)
{
	draw_map(dens_max);
}

static void draw_max_vel(void)
{
	draw_map(vel_max);
} */

static void draw_grid(void) 
{
	int i, j, i0, j0;
	float x, y, h, d00, d01, d10, d11;
	int gsize, slsize;
	float ratio = 0.3f; //0.05f;

	h = 1.0f/N;
	gsize = G;
	slsize = (N+2)/gsize;



	glBegin ( GL_QUADS );
		
		for (i = 0 ; i < G; i++) 
		for (j = 0 ; j < G; j++) {
			
			for (i0 = i * slsize; i0 <= (i+1)*slsize; i0++) {
				x = (i0-0.5f)*h;
				for (j0 = j * slsize; j0 <= (j+1)*slsize; j0++) {
					y = (j0-0.5f)*h;
					
					d00 = d01 = d10 = d11 = ratio * grid[i][j];

					glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
					glColor3f ( d10, d10, d10 ); glVertex2f ( x+h, y );
					glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
					glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+h );
				}
			}
		}
	
	glEnd ();
}


/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/

static void get_from_UI ( float ** d, float ** u, float ** v, int dir)
{
	int i, j, size = (N+2);

	for ( i=0 ; i<size ; i++ ) for( j=0 ; j<size ; j++) {
		u[i][j] = v[i][j] = d[i][j] = 0.0f;
		/* v[i] = force;  */ 
	}
/*
 * original way of obtaining forces and velocities
	if ( !mouse_down[0] && !mouse_down[2] ) return;

	i = (int)((       mx /(float)win_x)*N+1);
	j = (int)(((win_y-my)/(float)win_y)*N+1);

	if ( i<1 || i>N || j<1 || j>N ) return;

	if ( mouse_down[0] ) {
		u[IX(i,j)] = force * (mx-omx);
		v[IX(i,j)] = force * (omy-my);
	}

	if ( mouse_down[2] ) {
		d[IX(i,j)] = source;
	}
*/

	v[N/2][N/4] = force * dir; 
	d[N/2][N/4] = source * dir;

	omx = mx;
	omy = my;

	return;
}

/*
  ----------------------------------------------------------------------
   GLUT callback routines
  ----------------------------------------------------------------------
*/


static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}


/* 
 -------------------------------
  coordinating and comparing the computations in the two simulations
  ------------------------------- 
*/  

double compare_densities(float **a, float **b) {
	int i, j, size = N+2; 
	double sum = 0.0f; 
	for (i = 0; i < size; i++) 
	for (j = 0; j < size; j++) {
		sum += fabs(a[i][j] - b[i][j]);
	}
	sum /= size * size; 
	return sum;
}


static void step() {
	static int toggle = 0, dir = 1;
	// settle ui stuff
	if (toggle == 8) { 
			dir = 1 - dir; 
			toggle = 0 ;
	}
	toggle++;
	get_from_UI ( dens_prev, u_prev, v_prev, dir);
	get_from_UI ( denso_prev, uo_prev, vo_prev, dir);
	
	if (iter%20 == 0 ) { 
		refresh_grid(dens);
		dens_step = get_dens_step();
		vel_step = get_vel_step();
		cflag++;
	}
	
	//optim
	(*vel_step)( N, u, v, u_prev, v_prev, visc, dt);
	(*dens_step)( N, dens, dens_prev, u, v, diff, dt);
	//orig 
	(*vel_step_orig)( N, uo, vo, uo_prev, vo_prev, visc, dt);
	(*dens_step_orig)( N, denso, denso_prev, uo, vo, diff, dt);
	
	glutSetWindow ( win_id );
	glutPostRedisplay ();
	printf("ITER: %d\n", iter);
	printf("COMPARISON: %f\n", compare_densities(dens, denso));
	iter++;
}

static void key_func ( unsigned char key, int x, int y )
{
	
	switch ( key )
	{
		case 'c':
		case 'C':
			clear_data ();
			break;

		case 'q':
		case 'Q':
			free_data ();
			exit ( 0 );
			break;

		case 'v':
		case 'V':
			dvel++;
			dvel %= 3; 
			break;
		case 'p':
		case 'P':
			pause = (1 - pause); 
			break;
		case 'f':
		case 'F': 
			step();
			break;
	}
}

static void idle_func()
{
	if (pause) 
		return;
	step();
}

static void display_func ( void )
{
	pre_display ();

		switch(dvel) { 
			case 0: draw_density (); 
							break;
			case 1: draw_density_orig(); 
							break;
			case 2: draw_grid();
							break;
		}
	post_display ();
}


/*
  ----------------------------------------------------------------------
   open_glut_window --- open a glut compatible window and set callbacks
  ----------------------------------------------------------------------
*/

static void open_glut_window ( void )
{
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( win_x, win_y );
	win_id = glutCreateWindow ( "Alias | wavefront" );

	glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();
	glClear ( GL_COLOR_BUFFER_BIT );
	glutSwapBuffers ();

	pre_display ();

	glutKeyboardFunc ( key_func );
	glutReshapeFunc ( reshape_func );
	glutIdleFunc ( idle_func );
	glutDisplayFunc ( display_func );
}


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/
int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );

	if ( argc != 1 && argc != 8 ) {
		fprintf ( stderr, "usage : %s N dt diff visc force source grid\n", argv[0] );
		fprintf ( stderr, "where:\n" );\
		fprintf ( stderr, "\t N      : grid resolution\n" );
		fprintf ( stderr, "\t dt     : time step\n" );
		fprintf ( stderr, "\t diff   : diffusion rate of the density\n" );
		fprintf ( stderr, "\t visc   : viscosity of the fluid\n" );
		fprintf ( stderr, "\t force  : scales the mouse movement that generate a force\n" );
		fprintf ( stderr, "\t source : amount of density that will be deposited\n" );
		fprintf ( stderr, "\t grid : the size spot grid" ); 
		exit ( 1 );
	}

	if ( argc == 1 ) {
		N = 64;
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		force = 5.0f;
		source = 100.0f;
		G = 5; 
		fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force=%g source=%g\n",
			N, dt, diff, visc, force, source);
	} else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
		G = atoi(argv[7]);
	}

	printf ( "\n\nHow to use this demo:\n\n" );
	printf ( "\t Add densities with the right mouse button\n" );
	printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
	printf ( "\t Toggle diferent displays with the 'v' key\n" );
	printf ( "\t Clear the simulation by pressing the 'c' key\n" );
	printf ( "\t Quit by pressing the 'q' key\n" );

	dvel = 0;
	pause = 0;
	iter = 0;
	slot_size	= N / G;

	if ( !allocate_data () ) exit ( 1 );
	clear_data ();
	// DOM SOLVER STUFF
	cflag = 0;
	alloc_data();
	//original steps, for comparing with the optimized
	dens_step_orig = get_dens_step_orig();
	vel_step_orig = get_vel_step_orig();
	
	win_x = 512;
	win_y = 512;

	open_glut_window ();

	glutMainLoop ();

	exit ( 0 );
}
