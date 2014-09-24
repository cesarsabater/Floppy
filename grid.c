#include <stdio.h>
#include <stdlib.h>

#define max(x,y) (((x) > (y)) ? (x) : (y))

// extern variables
extern int N;
// global variables
int G;  // grid size 
int slot_size; 
int **grid, **grid_aux, **code_grid, **grid_new;	// grid of densities 

void alloc_matrix(int ***g)
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

void alloc_data() { 
	alloc_matrix(&grid);
	alloc_matrix(&grid_aux);
	alloc_matrix(&code_grid);
	alloc_matrix(&grid_new);
}

void grid_init()
{
	slot_size	= N / G;
	alloc_data();
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
	int orig_iter = N * N * 20;
 	int i, j, citer = 0; 
	for (i = 0; i < G; i++) 
	for (j = 0; j < G; j++) {
		citer += iter_from_level(code_grid[i][j]);
	}
	citer *= slot_size * slot_size;
	printf("executed iterations every lin_solve: %f (%d / %d)\n", (1.0 * citer) / orig_iter, citer, orig_iter);
}

int level_from_density(float dens) { 	
	if (dens > 0.5f) return 2;
	if (dens > 0.001f) return 1;
	return 0;  
}

int gridcmp(int **g1, int **g2) {
	int i, j, men = 0, size = G; 
	for (i = 0; i < size; i++) 
	for (j = 0; j < size; j++) {
		if (g1[i][j] != g2[i][j]) {  
			if (g1[i][j] < g2[i][j])
				men = -1;
			else 
				return 1;
		}
	}
	return men;
}

void gridcpy(int **src, int **dest) {
	int i, j, size = G; 
	for (i = 0; i < size; i++) 
	for (j = 0; j < size; j++) 
		dest[i][j] = src[i][j];
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
}
