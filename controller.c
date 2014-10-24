#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// simulation



/* global variables */
// simulation
int N = 200;
float dt=0.1, diff=0.00001, visc=0;
float **u, **v, **u_prev, **v_prev;
float **dens, **dens_prev;
extern void dens_step( int N, float **x, float **x0, float **u, float **v, float diff, float dt);
extern void vel_step( int N, float **u, float **v, float **u0, float **v0, float visc, float dt);
double time1, time2; 

//grid
int G = 10;  // grid size 
int slot_size; 
int **grid;

void alloc_grid(int ***g)
{
	int i;
	*g = (int**)malloc(G * sizeof(int*));
	if (!(*g)) 
		exit(0);
	for (i = 0; i < G; i++) {
		(*g)[i] = (int*)malloc(G * sizeof(int));
		if (!((*g)[i])) 
			exit(0);
	}
}	

void grid_init()
{
	slot_size	= N / G;
	alloc_grid(&grid);
}

void read_grid() { 
	int i, j;
	for (i = 0; i < G; i++) { 
		for (j = 0; j < G; j++) { 
			scanf("%d", &(grid[i][j]));
		}
	}
	for (i = 0; i < G; i++) { 
		for (j = 0; j < G; j++) { 
			printf("%d ", grid[i][j]);
		}
		printf("\n");
	}
}

float **alloc_matrix() { 
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

int allocate_data(void)
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

void read_matrix(float **m) {
	int i, j;
	int size = N+2;
	for(i = 0; i < size ; i++) {
		for(j = 0; j < size ; j++) {
			scanf("%f", &(m[i][j]));
		}
	}
}

void read_state() {
	read_matrix(dens);
	read_matrix(u);
	read_matrix(v);
}


int iter_from_level(int lev) { 
	switch (lev) { 
		case 2: 
			return 4; 
			break;
		case 1: 
			return 3; 
			break;
		case 0: 
			return 1;
			break;
		default: 
			return 20;
			break;
	}
}
	
int main() 
{
	grid_init();
	allocate_data();
	read_grid();
	read_state();
	dens_step( N, dens, dens_prev, u, v, diff, dt);
	vel_step( N, u, v, u_prev, v_prev, visc, dt);
	return 0;
}


