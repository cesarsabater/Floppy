#include <stdio.h>
#include <stdlib.h>

#define max(x,y) (((x) > (y)) ? (x) : (y))

// extern variables
extern int N;
// global variables
int G;  // grid size 
int slot_size; 
int **grid, **grid_aux, **code_grid;	// grid of densities 

void alloc_data() { 
	int i; 
	//alloc grid
	grid = (int**)malloc(G * sizeof(int*));
	grid_aux = (int**)malloc(G * sizeof(int*));
	code_grid = (int**)malloc(G * sizeof(int*));
	if (!grid || !grid_aux || !code_grid) 
		exit(0);
	for (i = 0; i < G; i++) {
		grid[i] = (int*)malloc(G * sizeof(int));
		grid_aux[i] = (int*)malloc(G * sizeof(int));
		code_grid[i] = (int*)malloc(G * sizeof(int));
		if (!(grid[i]) || !(grid_aux[i]) || !(code_grid[i])) 
			exit(0);
	}
}

void grid_init()
{
	slot_size	= N / G;
	alloc_data();
	//ensure that grids are different at initialize time
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

void calculate_iter() 
{
	int orig_iter = 200 * 200 * 20;
 	int i, j, citer = 0; 
	for (i = 0; i < G; i++) 
	for (j = 0; j < G; j++) {
		citer += iter_from_level(grid[i][j]);
	}
	citer *= slot_size * slot_size;
	printf("executed iterations every lin_solve: %f (%d / %d)\n", (1.0 * citer) / orig_iter, citer, orig_iter);
}

int level_from_density(float dens) { 	
	if (dens > 0.5f) return 2;
	if (dens > 0.001f) return 1;
	return 0;  
}

int compare_grids(int **g1, int **g2) {
	int i, j, size = G; 
	for (i = 0; i < size; i++) 
	for (j = 0; j < size; j++) {
		if (g1[i][j] != g2[i][j]) 
			return 0;
	}
	return 1;
}

/* 
 * TODO: check that this is OK!!!!
 */ 
void refresh_grid(float **dens) { 
	int i, j, i0, j0; 
	int gsize = G;  
	int size = N+2;
	int slotsize = size / gsize;
	int **aux;

	for (i = 0; i < gsize; i++) for (j = 0; j < gsize; j++) {
		grid[i][j] = 0;  
		// iterate in inside a slot
		for (i0 = i * slotsize; i0 < (i+1)*slotsize; i0++) 
		for (j0 = j * slotsize; j0 < (j+1)*slotsize; j0++)
			grid[i][j] = max(grid[i][j], level_from_density(dens[i0][j0]));
		// grid[i][j] /= 1.0f * (slotsize * slotsize); 
		//printf("grid[%d][%d]=%f\n", i, j, grid[i][j]);
	}
	for (i = 0; i < gsize; i++) for (j = 0; j < gsize; j++) {
		grid_aux[i][j] = grid[i][j]; 
		// 1st comparison
		if (i > 0 && j > 0) { 
			if 			(grid_aux[i-1][j-1] < grid[i][j]) grid_aux[i-1][j-1] = grid[i][j];
			else if (grid_aux[i][j] < grid[i-1][j-1]) grid_aux[i][j] = grid[i-1][j-1];
		}
		// 2nd comparison
		if (i > 0) { 
			if   		(grid_aux[i-1][j] < grid[i][j]) 	grid_aux[i-1][j] = grid[i][j]; 
			else if (grid_aux[i][j] < grid[i-1][j])  	grid_aux[i][j] = grid[i-1][j];
		}
		// 3rd comparison
		if (j > 0) { 
			if 			(grid_aux[i][j-1] < grid[i][j]) 	grid_aux[i][j-1] = grid[i][j];
			else if (grid_aux[i][j] < grid[i][j-1]) 	grid_aux[i][j] = grid[i][j-1];
		}
		// 4th comparison
		if (j > 0 && i < gsize-1) {
			if 			(grid_aux[i+1][j-1] < grid[i][j]) grid_aux[i+1][j-1] = grid[i][j];
			else if (grid_aux[i][j] < grid[i+1][j-1]) grid_aux[i][j] = grid[i+1][j-1];
		}
	} 
	aux = grid;
	grid = grid_aux;
	grid_aux = aux;
	calculate_iter();
}
