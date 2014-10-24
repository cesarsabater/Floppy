#include <stdio.h>
#include <stdlib.h>

#define max(x,y) (((x) > (y)) ? (x) : (y))

// extern variables
extern int N;
// global variables
int G;  // grid size 
int slot_size; 
int **grid, **grid_aux;
int *bounds;

/////////////////////////////////////////////////
/////////// USER PROVIDED INFORMATION ///////////
/////////////////////////////////////////////////

// que iteradores son el espacio de iteracion??

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

int level_from_density(float dens) { 	
	if (dens > 0.5f) return 2;
	if (dens > 0.001f) return 1;
	return 0;  
}

int get_levels() 
{
	return 3;
}

int max_outer_loop_iterations() {
	return 20;
}

/////////////////////////////////////////////////
/////// END OF USER PROVIDED INFORMATION ////////
/////////////////////////////////////////////////

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

static void free_matrix(int **m) {
	int i;
	for (i = 0; i < G; i++) 
		if (m[i]) 
			free(m[i]);
	free(m);
	if (bounds)
		free(bounds);
}

void alloc_data() { 
	alloc_matrix(&grid);
	alloc_matrix(&grid_aux);
	bounds = (int*)malloc(G * sizeof(int));
}

void free_data() {
	free_matrix(grid);
	free_matrix(grid_aux);
}

void finish_grid() {
	free_data();
}

int lower_bound(i) { 
	return i*slot_size + 1;
}

void init_bounds() {
	int i;
	for (i = 0; i < G; i++) { 
		bounds[i] = i*slot_size + 1;
	}
}

int get_slot(int i) { 
	return (i-1)/slot_size;
}

void grid_init()
{
	G = 10;
	slot_size	= N / G;
	alloc_data();
	//init_bounds();
}

void calculate_iter(int **grid) 
{
	int orig_iter = N * N * 20;
 	int i, j, citer = 0; 
	for (i = 0; i < G; i++) 
	for (j = 0; j < G; j++) {
		citer += iter_from_level(grid[i][j]);
	}
	citer *= slot_size * slot_size;
	printf("executed iterations every lin_solve: %f (%d / %d)\n", (1.0 * citer) / orig_iter, citer, orig_iter);
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

