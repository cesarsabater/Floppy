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
#include <omp.h> 

/* macros */

#define DISPLAY_DUMPFILE "fluid_orig_dump.out"
#define DISPLAY_DUMPRATE 100

#define IX(i,j) ((i)+(N+2)*(j))
#define max(x,y) (((x) > (y)) ? (x) : (y))



// simulation
extern void vel_step ( int N, float **u, float **v, float **u0, float **v0, float visc, float dt);
void dens_step ( int N, float **x, float **x0, float **u, float **v, float diff, float dt);
/* global variables */
int N, NUM_ITER;
int pause;
int iter;
static float dt, diff, visc;
static float force, source;
static int dvel;

float ** u, ** v, ** u_prev, ** v_prev;
float ** dens, ** dens_prev;
float ** vel_max, ** dens_max;

static int win_id;
static int win_x, win_y;
static int omx, omy, mx, my;

FILE *dumpfile;
double time1, time2; 
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
}

static void clear_data (void)
{
	int i, j, size=N+2;

	for ( i=0 ; i<size ; i++ ) for (j=0 ; j<size; j++) {
		dens_max[i][j] = vel_max[i][j] = u[i][j] = v[i][j] = u_prev[i][j] = v_prev[i][j] = dens[i][j] = dens_prev[i][j] = 0.0f;
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

	u				  = alloc_matrix(); v 			  = alloc_matrix();
	u_prev	  = alloc_matrix(); v_prev	  = alloc_matrix();
	vel_max   = alloc_matrix();
	dens		  = alloc_matrix(); dens_prev	= alloc_matrix();
	dens_max  = alloc_matrix();
	
	if ( !u || !v || !u_prev || !v_prev ||!vel_max || !dens || !dens_prev || !dens_max) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}
	
	return ( 1 );
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

static void draw_density ( void )
{
	draw_map(dens);
}

static void draw_max_dens(void)
{
	draw_map(dens_max);
}

static void draw_max_vel(void)
{
	draw_map(vel_max);
}


void dump_matrix(float **m, FILE *df) {
	int i, j;
	fprintf(df, "%d\n", iter); 
	for(i = 0; i <= N ; i++) {
		for(j = 0; j <= N ; j++) 
			fprintf(df, "%f ", m[i][j]);
		fprintf(df, "\n");
	}
}
/*
  ----------------------------------------------------------------------
   relates mouse movements to forces sources
  ----------------------------------------------------------------------
*/
static void get_from_UI ( float ** d, float ** u, float ** v )
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
	static int toggle = 0, dir = 1;
	if (toggle == 8) { dir = ( dir == 0 ) ? 1 : 0; toggle = 0 ;}
	toggle++;
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

static void poll_dens(int N, float **d, float **dm)  
{
	int i,j;
	for ( i=0 ; i<=N ; i++ ) 
	for ( j=0 ; j<=N ; j++ ) {
		dm[i][j] = max(dm[i][j], d[i][j]);
	}
}

static void poll_vel(int N, float **u, float **v, float **vm)
{
	int i,j;
	for ( i=0 ; i<=N ; i++ ) 
	for ( j=0 ; j<=N ; j++ ) {
		vm[i][j] = max(vm[i][j], sqrt((u[i][j]*u[i][j]) +  (v[i][j]*v[i][j])) * 10);
	}
}


void finish_sim() {
	free_data();
	fclose(dumpfile);
}

static void step() {
	get_from_UI ( dens_prev, u_prev, v_prev );
	
	vel_step(N, u, v, u_prev, v_prev, visc, dt);
	dens_step(N, dens, dens_prev, u, v, diff, dt);
	//poll_dens(N, dens, dens_max); 
	//poll_vel(N, u, v, vel_max); 
	
	glutSetWindow ( win_id );
	glutPostRedisplay ();
	//printf("ITER: %d\n", iter);
	iter++;
	if (iter >= NUM_ITER) {
		time2 = omp_get_wtime();
    printf( "Number of seconds: %f\n", time2 - time1);
		finish_sim();
		exit(0);
	}
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
			finish_sim();
			exit ( 0 );
			break;
		case 'v':
		case 'V':
			dvel++;
			dvel %= 4; 
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

static void idle_func ( void )
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
			case 1: draw_velocity (); 
							break;
			case 2: draw_max_dens (); 
							break;
			case 3: draw_max_vel (); 
							break;
			
	}
	
	if (iter % DISPLAY_DUMPRATE == 0) 
		dump_matrix(dens, dumpfile); 
	
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
	//glutMouseFunc ( mouse_func );
	//glutMotionFunc ( motion_func );
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
		fprintf ( stderr, "usage : %s N dt diff visc force source\n", argv[0] );
		fprintf ( stderr, "where:\n" );\
		fprintf ( stderr, "\t N      : grid resolution\n" );
		fprintf ( stderr, "\t dt     : time step\n" );
		fprintf ( stderr, "\t diff   : diffusion rate of the density\n" );
		fprintf ( stderr, "\t visc   : viscosity of the fluid\n" );
		fprintf ( stderr, "\t force  : scales the mouse movement that generate a force\n" );
		fprintf ( stderr, "\t source : amount of density that will be deposited\n" );
		fprintf ( stderr, "\t num_iter : number of iterations to be performed\n" );
		exit ( 1 );
	}

	if ( argc == 1 ) {
		N = 64;
		dt = 0.1f;
		diff = 0.0f;
		visc = 0.0f;
		force = 5.0f;
		source = 100.0f;
		NUM_ITER = 600;
		fprintf ( stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force=%g source=%g num_iter=%d\n",
			N, dt, diff, visc, force, source, NUM_ITER);
	} else {
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		diff = atof(argv[3]);
		visc = atof(argv[4]);
		force = atof(argv[5]);
		source = atof(argv[6]);
		NUM_ITER = atoi(argv[7]);
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

	if ( !allocate_data () ) exit ( 1 );
	clear_data ();
	
	// DOM SOLVER STUFF
	dumpfile = fopen(DISPLAY_DUMPFILE, "w");


	win_x = 512;
	win_y = 512;

	open_glut_window ();

	time1 = omp_get_wtime();
	
	glutMainLoop ();

	exit ( 0 );
}
