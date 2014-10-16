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
#include <pthread.h>
#include <omp.h>

/* macros */

#define DISPLAY_DUMPFILE "fluid_dump.out"
#define DISPLAY_DUMPRATE 100
#define max(x,y) (((x) > (y)) ? (x) : (y))

/* external definitions  */
//generator
extern void transformer_init();
extern void get_new_code_if_available();
extern void start_generator(void*);
//grid
extern int G;
extern int **grid, **code_grid, **grid_aux;
extern int slot_size;
extern void grid_init();
extern void refresh_grid(float**);
extern int gridcmp(int**, int**);
extern void calculate_iter();
// main
extern int N, NUM_ITER;
extern float dt, diff, visc;
extern float force, source;
extern pthread_mutex_t fmutex;
extern pthread_mutex_t gmutex;
// original simulation
extern void dens_step (int, float**, float**, float**, float**, float, float);
extern void vel_step (int, float**, float**, float**, float**, float, float);

/* global variables */

static int dvel;
int pause;
int iter;

float ** u, ** v, ** u_prev, ** v_prev;
float ** dens, ** dens_prev;

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

static void free_data(void)
{
	if(u) free_matrix(u);
	if(v) free_matrix(v);
	if(u_prev) free_matrix(u_prev);
	if(v_prev) free_matrix(v_prev);
	if(dens) free_matrix(dens);
	if(dens_prev) free_matrix(dens_prev);
}

static void clear_data (void)
{
	int i, j, size=N+2;

	for ( i=0 ; i<size ; i++ ) for (j=0 ; j<size; j++) {
		u[i][j] = v[i][j] = u_prev[i][j] = v_prev[i][j] = dens[i][j] = dens_prev[i][j] = 0.0f;
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
	dens		  = alloc_matrix(); dens_prev	= alloc_matrix();
	if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
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

void dump_matrix(float **m) {
	int i, j;
	fprintf(dumpfile, "%d\n", iter); 
	for(i = 0; i <= N ; i++) {
		for(j = 0; j <= N ; j++) 
			fprintf(dumpfile, "%f ", m[i][j]);
		fprintf(dumpfile, "\n");
	}
}

static void draw_density(void)
{
	draw_map(dens);
}

static void draw_grid(void) 
{
	// TODO: please optimize this mess!
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
					
					d00 = d01 = d10 = d11 = ratio * code_grid[i][j];

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

void apply_sim() 
{
	//mutex lock
	/*
	pthread_mutex_lock(&fmutex);
	if (dens_step_opt && vel_step_opt && gridcmp(grid_aux, code_grid) < 1) {
		(*vel_step_opt)( N, u, v, u_prev, v_prev, visc, dt);
		(*dens_step_opt)( N, dens, dens_prev, u, v, diff, dt);
		//mutex unlock and exit
		calculate_iter();
		pthread_mutex_unlock(&fmutex); 
		return;
	}
	//mutex unlock
	pthread_mutex_unlock(&fmutex);*/
	
	pthread_mutex_lock(&fmutex);
	vel_step( N, u, v, u_prev, v_prev, visc, dt);
	dens_step( N, dens, dens_prev, u, v, diff, dt);
	if (iter % DISPLAY_DUMPRATE == 0) 
		calculate_iter(); 
	pthread_mutex_unlock(&fmutex);
	//printf("Run original code: NO ITERATIONS SAVED!\n");
}

void close_display() { 
	free_data();
	fclose(dumpfile);
}

static void step() {
	//refresh grid
	//pthread_mutex_lock(&gmutex);
	refresh_grid(dens);
	//pthread_mutex_unlock(&gmutex);
	//get density and velocity new values
	get_from_UI(dens_prev, u_prev, v_prev);
	//apply simulation step
	apply_sim();
	//display
	glutSetWindow ( win_id );
	glutPostRedisplay ();
	//printf("ITER: %d\n", iter);
	iter++;
	if (iter >= NUM_ITER) {
		time2 = omp_get_wtime();
    printf( "Number of seconds: %f\n", time2 - time1 );
		close_display();
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
			close_display(); 
			exit(0);
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

static void display_func(void)
{
	pre_display();
	switch(dvel) { 
		case 0: draw_density(); 
						break;
		case 1: draw_velocity(); 
						break;
		case 2: draw_grid();
						break;
	}
	
	if (iter % DISPLAY_DUMPRATE == 0) 
		dump_matrix(dens);
	
	post_display();
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

void display_init(int argc, char **argv) {
	/*
	dens_step_opt = 0;
	vel_step_opt = 0;
	*/
	glutInit(&argc, argv);
	dvel = 0;
	pause = 0;
	iter = 0;
	if (!allocate_data()) 
		exit (1);
	clear_data();
	win_x = 512;
	win_y = 512;
	
	dumpfile = fopen(DISPLAY_DUMPFILE, "w");
	
	open_glut_window();
}
/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/
void start_sim(void *arg)
{
	pthread_t generator;
	//one step to put values to the grid
	
	step();
	//grid is ready to be used by the code generator
	pthread_create(&generator, NULL, (void*)start_generator, NULL);
	//start main simulation loop
	time1 = omp_get_wtime();
	glutMainLoop();
	pthread_join(generator, NULL);
}
