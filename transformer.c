// Compile with -ldl
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dlfcn.h>


#include <isl/ctx.h>
#include <isl/options.h>
#include <isl/printer.h>
#include <isl/set.h>

#include <spot/spot.h>
#include <osl/osl.h>
#include <osl/extensions/doi.h>

#define DENS_FUNCTION "dens_step"
#define VEL_FUNCTION "vel_step"
#define FILENAME     "solver3.c"
#define BASE_CODE "simulation_base.c"
#define ORIGINAL_CODE "simulation_original.c"
#define GENERATED_CODE "simulation_generated.c"
#define COMPILED 		 "libsimulator.so"
#define COMPILED_ORIG "libsimulator_orig.so"
#define LEVELS 3
#define MAX_ITERATIONS 20
#define FILENAME_LEN 20d

#define max(x,y) (((x) > (y)) ? (x) : (y))

/* extern variables */
extern int N;

extern float ** u, ** v, ** u_prev, ** v_prev;
extern float ** dens, ** dens_prev;
extern float ** vel_max, ** dens_max;

/* local variables */ 
int G;  /* grid size */
int slot_size; 

int **grid, **grid_aux;	/* grid of densities */ 
osl_scop_p scop; 
osl_spot_p spots; 
int cflag;

isl_ctx *ctx;
isl_set **sets; 

void *handle = NULL, *handle_orig = NULL;

void gen_random(char *s, const int len) {
    int i;
    static const char alphanum[] =
        "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for (i = 0; i < len; ++i) 
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    s[len] = '\0';
}

void alloc_data() { 
	int i; 
	//alloc grid
	grid = (int**)malloc(G * sizeof(int*));
	grid_aux = (int**)malloc(G * sizeof(int*));
	if (!grid || !grid_aux) 
		exit(0);
	for (i = 0; i < G; i++) {
		grid[i] = (int*)malloc(G * sizeof(int));
		grid_aux[i] = (int*)malloc(G * sizeof(int));
		if (!(grid[i]) || !(grid_aux[i])) 
			exit(0);
	}
}

void create_isl_set_head(char **params, int nparams, char **iter, int niter, 
															char *s, int *high_water_mark) { 
	int i;
	//open params 
	s[0] = '\0';
	if (nparams > 0 && params) { 
		osl_util_safe_strcat(&s, "[", high_water_mark); 
		osl_util_safe_strcat(&s, params[0], high_water_mark);
		for (i=1; i < nparams; i++) { 
			osl_util_safe_strcat(&s, ",", high_water_mark);
			osl_util_safe_strcat(&s, params[i], high_water_mark);
		}
		osl_util_safe_strcat(&s, "]->{", high_water_mark); 
	}
	if (niter > 0 && iter) { 
		osl_util_safe_strcat(&s, "[", high_water_mark);  
		osl_util_safe_strcat(&s, iter[0], high_water_mark);  
		for (i=1; i < niter; i++) { 
			osl_util_safe_strcat(&s, ",", high_water_mark);
			osl_util_safe_strcat(&s, iter[i], high_water_mark);
		} 
		osl_util_safe_strcat(&s, "]:", high_water_mark);  
	} 
}

int level_from_density(float dens) { 	
	if (dens > 0.5f) return 2;
	// return 2; 
	//return 3;
	if (dens > 0.001f) return 1;
	//return 1;
	//return 0;  
	return 0;  
}

int iter_from_level(int lev) { 
	switch (lev) { 
	/*	case 3:
			return 4;
			break;*/
		case 2: 
			return 4; 
			break;
		case 1: 
			return 3; 
			break;
		case 0: 
			return 1;
			break;
	}
	return 0;
	//return 2;
	//return 3;
	int calc = ((lev+1)*MAX_ITERATIONS) / (LEVELS);
/*	if (calc <= 0) 
		calc = 2;  */
	return calc;
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

void create_isl_set_from_grid() {
	char **scop_params, **scop_iter;
	int niter, i, j; 
	osl_body_p stmt_body;
	char buffer[OSL_MAX_STRING];
	int high_water_mark = OSL_MAX_STRING;
	char *s;
	
	ctx = isl_ctx_alloc();		 
	assert(ctx);
	isl_options_set_on_error(ctx, ISL_ON_ERROR_ABORT);
	
	// get parameters
	if(scop->context->nb_parameters) 
		scop_params = ((osl_strings_p)scop->parameters->data)->string;
	// get iterators
	niter = osl_statement_get_nb_iterators(scop->statement);
	stmt_body = osl_statement_get_body(scop->statement);
	if(scop->context->nb_parameters)
		scop_iter = stmt_body->iterators->string;
	
	
	// alloc and init data
	s = (char*)malloc(high_water_mark * sizeof(char)); 
	sets = (isl_set**)malloc(LEVELS * sizeof(isl_set*));
	for (i = 0; i < LEVELS; i++) 
		sets[i] = NULL;
	
	
	// compute domains
	for (i = 0; i < G; i++)  
	for (j = 0; j < G; j++) {
			int level = grid[i][j];
			isl_set *set_i;
			s[0] ='\0';
			high_water_mark = OSL_MAX_STRING;
			// construct isl set string
			create_isl_set_head(scop_params, scop->context->nb_parameters, 
														scop_iter, niter, s, &high_water_mark);														
			// grid domain constraints
			sprintf(buffer, "%s > %d", scop_iter[1], (int)(i*slot_size)); 
			osl_util_safe_strcat(&s, buffer, &high_water_mark);
			sprintf(buffer, " and %s <= %d", scop_iter[1], (int)((i+1)*slot_size)); 
			osl_util_safe_strcat(&s, buffer, &high_water_mark);
			sprintf(buffer, " and %s > %d", scop_iter[2], (int)(j*slot_size)); 
			osl_util_safe_strcat(&s, buffer, &high_water_mark);
			sprintf(buffer, " and %s <= %d}", scop_iter[2], 1+(int)((j+1)*slot_size)); 
			osl_util_safe_strcat(&s, buffer, &high_water_mark);
			
			// read domain
			set_i = isl_set_read_from_str(ctx, s);
			if (!sets[level])
				sets[level] = isl_set_copy(set_i);
			else 
				sets[level] = isl_set_union(sets[level], set_i);
	}
	// add iteration constraints and simplify sets 
	for (i = 0; i < LEVELS; i++) { 
		if (sets[i]) {
			int iter = iter_from_level(i);
			isl_set *set_k;
			s[0] = '\0';
			high_water_mark = OSL_MAX_STRING;
			// construct isl set string
			create_isl_set_head(scop_params, scop->context->nb_parameters, 
														scop_iter, niter, s, &high_water_mark);	
			// iteration constraints
			sprintf(buffer, "%s >= %d and %s < %d}", scop_iter[0], iter, scop_iter[0], MAX_ITERATIONS);
			osl_util_safe_strcat(&s, buffer, &high_water_mark);
			// read domain
			set_k = isl_set_read_from_str(ctx, s);
			sets[i] = isl_set_intersect(sets[i], set_k);
			sets[i] = isl_set_coalesce(sets[i]);
			
			printf("LEVEL %d\n", i);
			isl_set_print(sets[i], stdout, 0, ISL_FORMAT_ISL);
			printf("\n");
		}
	}
	// free stuff
	free(s);
	// TODO: FREE CONTEXT
}

void create_spots() 
{ 
	int i;
	char *str;
	isl_printer *p = NULL;
	osl_spot_p spot;
	osl_generic_p generic;
	
	p = isl_printer_to_str(ctx);
	p = isl_printer_set_output_format(p, ISL_FORMAT_ISL); 
	// remove old spots
	osl_generic_remove(&(scop->extension), OSL_URI_SPOT);
	// add new spots
	for (i = 0; i < LEVELS; i++)
		if (sets[i]) { 
			spot = osl_spot_malloc();
			spot->priority = i; 
			spot->comp = NULL; //osl_util_strdup("\""OSL_SPOT_NULL"\"");
			//spot->user = isl_set_copy(sets[i]);
			p = isl_printer_print_set(p, sets[i]);
			str = isl_printer_get_str(p);
			p = isl_printer_flush(p); 
			spot->dom = osl_util_strdup(str);
			spots = osl_spot_concat(spot, spots);
		}
		
		
	if (spots) { 
		generic = osl_generic_shell(spots, osl_spot_interface());
		osl_generic_add(&scop->extension, generic);
	}

	// free sets??

	// free printer 
	isl_printer_free(p);	
}

void gen_code() {
	FILE *input, *output; 
	printf("generacion de codigo!\n");
	// read file	
	input = fopen(BASE_CODE, "r");
	scop = spot_scop_read_from_c(input, BASE_CODE);
	fclose(input);
	if (!scop) {
		fprintf(stderr, "[dom_solver.c] NULL SCOP (test)\n");
		exit(1); 
	}
	create_isl_set_from_grid();  
	create_spots();
	spot_compute_scops(scop);
	
	printf("\n\n\nSPOT\n");
	//osl_scop_dump(stdout, scop);
	
	output = fopen(GENERATED_CODE, "w");
	spot_scop_print_to_c(output, scop);
	 /*if (output != NULL) 
		fclose(output);*/
	// DEBUG
	
	//osl_spot_dump(stdout, spots);
	
	//osl_spot_free(spots);
	spots = NULL;
	osl_scop_free(scop);
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
	
	gen_code();
}

/**
 * 
 * PASOS A SEGUIR: 
 * 
 * 1- Generar grilla de computaciones
 * 
 * 2- Genearar dominios isl a partir de la grilla
 * 		a- crear funciones spot para manejar la libreria programaticamente
 * 		b- computar los nuevos dominios
 * 		c- insertarlos en doi
 *
 * 3- Recompilar el codigo y hacerlo correr
 * 
 * 4- Hacer correr el thread que ahce esto en paralelo a la simulacion
 * 
 */ 

typedef void (*stepfun)(int,float**,float**,float**,float**,float,float);
/*
stepfun compile_and_load_step(void **handle, char *source, char *fun, char *obj) { 
	stepfun f;
  char *error;
  char command[50]; 
  printf("compilacion!\n");
  
  if (!(handle)) { 
		
		// compile the code
		sprintf(command, "gcc -fPIC -shared -o %s %s", obj, source);
		system(command);
		// load the library
		
		handle = dlopen(obj, RTLD_LAZY);
		if (!(handle)) {
			fprintf(stderr, "%s\n", dlerror());
			exit(EXIT_FAILURE);
		}
	}
	
	// Get the pointer to the function we want to execute
	f = (stepfun)dlsym(handle, fun);
	// check for errors
	if ((error = dlerror()) != NULL) {
    fprintf(stderr, "%s\n", error);
    exit(EXIT_FAILURE);
  }
  
  printf("111\n");
	return f;
} */

/*
stepfun get_step(int orig, char *fun) {
	stepfun f;
	char buffer[50];
	 char *error;
	//printf("CFLAG: %d\n", cflag);
	if (handle != NULL) 
			dlclose(handle); 
	if (cflag && orig == 0) {
		sprintf(buffer, "sed -i 's/versionXXXX/%d/g' "GENERATED_CODE, cflag);
		printf("%s\n", buffer);
		system(buffer);
		sprintf(buffer, "gcc -fPIC -shared -o %s %s", COMPILED, GENERATED_CODE);
		system(buffer);
		handle = dlopen(COMPILED, RTLD_LAZY);
		f = (stepfun)dlsym(handle, fun);
		//f = compile_and_load_step(&handle, GENERATED_CODE, fun, COMPILED);
	}
	else { 
		system("cp "BASE_CODE" "ORIGINAL_CODE); 
		sprintf(buffer, "sed -i 's/versionXXXX/%d/g' "ORIGINAL_CODE, 0);
		printf("%s\n", buffer);
		system(buffer);
		sprintf(buffer, "gcc -fPIC -shared -o %s %s", COMPILED_ORIG, ORIGINAL_CODE);
		system(buffer);
		handle = dlopen(COMPILED_ORIG, RTLD_LAZY);
		//f = (stepfun)dlsym(handle_orig, fun);
		
		//f = compile_and_load_step(&handle_orig, ORIGINAL_CODE, fun, COMPILED_ORIG);
		f = (stepfun)dlsym(handle, fun);
	}	
	// Get the pointer to the function we want to execute
	//f = (stepfun)dlsym(handle, fun);
	// check for errors
	if ((error = dlerror()) != NULL) {
    fprintf(stderr, "%s\n", error);
    exit(EXIT_FAILURE);
  }
	
	return f;
}*/


stepfun get_step(int orig, char *fun) {
  void (*f)(int,float**,float**,float**,float**,float,float);
  char buffer[50]; 
  char *error;
  printf("compilacion!\n");
  if (handle != NULL) 
		dlclose(handle); 
	// dynamic compilation of some code (which can be dynamically generated!)
	printf("CFLAG: %d\n", cflag);
  if (cflag && orig == 0) {
		sprintf(buffer, "sed -i 's/versionXXXX/%d/g' "GENERATED_CODE, cflag);
		printf("%s\n", buffer);
		system(buffer);
		system("gcc -fPIC -shared -o "COMPILED" "GENERATED_CODE);
		printf("gcc -fPIC -shared -o "COMPILED" "GENERATED_CODE"\n");
		handle = dlopen(COMPILED, RTLD_LAZY);
		// load the dynamic library which have been dynamically generated
		if (!handle) {
			fprintf(stderr, "%s\n", dlerror());
			exit(EXIT_FAILURE);
		}
		// Get the pointer to the function we want to execute
		f = (void (*)(int,float**,float**,float**,
											float**,float,float))dlsym(handle, fun);
	}
	else { 
		system("cp "BASE_CODE" "ORIGINAL_CODE); 
		sprintf(buffer, "sed -i 's/versionXXXX/%d/g' "ORIGINAL_CODE, cflag);
		system(buffer);
		system("gcc -fPIC -shared -o "COMPILED_ORIG" "ORIGINAL_CODE);
		 handle_orig = dlopen(COMPILED_ORIG, RTLD_LAZY);
		 // load the dynamic library which have been dynamically generated
		if (!handle_orig) {
			fprintf(stderr, "%s\n", dlerror());
			exit(EXIT_FAILURE);
		}
		// Get the pointer to the function we want to execute
		f = (void (*)(int,float**,float**,float**,
										float**,float,float))dlsym(handle_orig, fun);
	}
	
  if ((error = dlerror()) != NULL) {
    fprintf(stderr, "%s\n", error);
    exit(EXIT_FAILURE);
  }
  
	return f;
}

stepfun get_dens_step_orig() {
	return get_step(1, DENS_FUNCTION); 
}

stepfun get_vel_step_orig() {
	return get_step(1, VEL_FUNCTION); 
}

stepfun get_dens_step() {
	return get_step(0, DENS_FUNCTION);
}

stepfun get_vel_step(){
	return get_step(0, VEL_FUNCTION);
}


