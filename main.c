#include <stdio.h>
#include <pthread.h>
#include <GL/glut.h>

/* extern */
//grid
extern int **grid, **grid_aux, **old_grid;
extern int G;
extern void grid_init();
//display
extern void display_init(int, char**);
extern void start_sim(void*);
//generator
extern void start_generator(void*);
extern void transformer_init();
/* global */
int N;
float dt, diff, visc;
float force, source;
pthread_mutex_t fmutex = PTHREAD_MUTEX_INITIALIZER;

void init(int argc, char ** argv) {
	display_init(argc, argv);
	transformer_init();
	grid_init();
}

int main (int argc, char ** argv)
{
	printf("1!\n");
	pthread_t generator, simulator;
	
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

	printf("2!\n");

	init(argc, argv);
	
	printf("3!\n");

	// exec threads here
	//pthread_create(&simulator, NULL, (void*)start_sim, NULL);
	pthread_create(&generator, NULL, (void*)start_generator, NULL);
	//printf("4!\n");
	//pthread_join(simulator, NULL);

	//printf("5!\n");
	start_sim(NULL);
	pthread_join(generator, NULL);
	
	exit ( 0 );
}
