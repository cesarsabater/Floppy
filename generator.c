// Compile with -ldl
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dlfcn.h>
#include <pthread.h>

#include <isl/ctx.h>
#include <isl/options.h>
#include <isl/printer.h>
#include <isl/set.h>

#include <spot/spot.h>
#include <osl/osl.h>
#include <osl/extensions/doi.h>

/*
#define DENS_FUNCTION 	"dens_step"
#define VEL_FUNCTION 		"vel_step"
*/
#define OPT_FUNCTION 		"lin_solve_opt"
#define BASE_CODE 			"lin_solve_base.c"
#define GENERATED_CODE 	"lin_solve_generated.c"
#define COMPILED 		 		"libsimulator.so"
#define COMPILED_ORIG 	"libsimulator_orig.so"
#define FILENAME_LEN 20

/* extern variables */
//main
extern int N;
extern pthread_mutex_t fmutex;
extern pthread_mutex_t gmutex;
//display
extern float ** u, ** v, ** u_prev, ** v_prev;
extern float ** dens, ** dens_prev;
/*
extern void (*vel_step_opt)(int, float**,float**, float**, float**, float, float);
extern void (*dens_step_opt)(int, float**, float**, float**, float**, float, float);
*/
//grid
extern int G;  // grid size 
extern int slot_size; 
extern int **grid, **grid_aux, **code_grid, **grid_new;	// grid of densities 
extern void refresh_grid(float**);
extern int gridcmp(int **, int **);
extern int gridcpy(int **, int **);
extern int iter_from_level(int);
extern int get_levels();
extern int max_outer_loop_iterations();
//simulation_original
extern void lin_solve_original(int, int, float **, float **, float, float);
/* local variables */ 
int cflag;
void *handle = NULL;
char obj_name[50];
osl_scop_p scop; 
osl_spot_p spots; 
isl_ctx *ctx;
isl_set **sets; 
void (*lin_solve_opt)(int, int, float **, float **, float, float);

//typedef void (*stepfun)(int,float**,float**,float**,float**,float,float);
typedef void (*stepfun)(int, int , float **, float **, float, float);

void lin_solve_safe( int N, int b, float **x, float **x0, float a, float c) { 
	// use de mutex here or before in sumulation-orig
	//pthread_mutex_lock(&fmutex);
	if (lin_solve_opt && gridcmp(grid_aux, code_grid) < 1) {
		(*lin_solve_opt)(N, b, x, x0, a, c);
		//mutex unlock and exit
		//calculate_iter();
		//pthread_mutex_unlock(&fmutex); 
		return;
	}
	//mutex unlock
	//pthread_mutex_unlock(&fmutex);
	lin_solve_original(N, b, x, x0, a, c);
	printf("Run original code: NO ITERATIONS SAVED!\n");
}

void gen_random_obj(char *s, int len) {
    int i;
    static const char alphanum[] =
        "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for (i = 0; i < len-3; ++i) 
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    strcpy(s+len-3, ".so");
    s[len] = '\0';
}

void transformer_init() { 
	lin_solve_opt = NULL;
	cflag = 0;
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
	sets = (isl_set**)malloc(get_levels() * sizeof(isl_set*));
	for (i = 0; i < get_levels(); i++) 
		sets[i] = NULL;
	// compute domains
	for (i = 0; i < G; i++)  
	for (j = 0; j < G; j++) {
			int level = grid_new[i][j];
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
	for (i = 0; i < get_levels(); i++) { 
		if (sets[i]) {
			int iter = iter_from_level(i);
			isl_set *set_k;
			s[0] = '\0';
			high_water_mark = OSL_MAX_STRING;
			// construct isl set string
			create_isl_set_head(scop_params, scop->context->nb_parameters, 
														scop_iter, niter, s, &high_water_mark);	
			// iteration constraints
			sprintf(buffer, "%s >= %d and %s < %d}", scop_iter[0], iter, scop_iter[0], max_outer_loop_iterations());
			osl_util_safe_strcat(&s, buffer, &high_water_mark);
			// read domain
			set_k = isl_set_read_from_str(ctx, s);
			sets[i] = isl_set_intersect(sets[i], set_k);
			sets[i] = isl_set_coalesce(sets[i]);
			//DEBUG
			/*
			printf("LEVEL %d\n", i);
			//isl_set_print(sets[i], stdout, 0, ISL_FORMAT_ISL);
			printf("\n");
			*/
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
	for (i = 0; i < get_levels(); i++) {
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
	//printf("generacion de codigo!\n");
	// read file	
	input = fopen(BASE_CODE, "r");
	scop = spot_scop_read_from_c(input, BASE_CODE);
	fclose(input);
	if (!scop) {
		fprintf(stderr, "[dom_solver.c] NULL SCOP (gen_code)\n");
		exit(1); 
	}
	create_isl_set_from_grid();  
	create_spots();
	spot_compute_scops(scop);
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

void get_steps(stepfun *l) {
  char buffer[50]; 
  char *error;
  //printf("compilacion!\n");
	// dynamic compilation of some code (which can be dynamically generated!)
	sprintf(buffer, "sed -i 's/versionXXXX/%d/g' "GENERATED_CODE, cflag);
	system(buffer);
	
	system("tcc -fPIC -shared -o "COMPILED" "GENERATED_CODE);
	
	// CONCURRENT REGION
	pthread_mutex_lock(&fmutex);
	if (handle != NULL) 
			dlclose(handle);  
	handle = dlopen(COMPILED, RTLD_LAZY);
	// TODO: CHECK THIS!
	// load the dynamic library which have been dynamically generated
	if (!handle) {
		fprintf(stderr, "%s\n", dlerror());
		exit(EXIT_FAILURE);
	}
	// Get the pointer to the function we want to execute
	*l = (void (*)(int, int , float **, float **, float, float))dlsym(handle, OPT_FUNCTION);
	/*
	*d = (void (*)(int,float**,float**,float**,
										float**,float,float))dlsym(handle, DENS_FUNCTION);
	*v = (void (*)(int,float**,float**,float**,
										float**,float,float))dlsym(handle, VEL_FUNCTION);
	*/ 
	// refresh code grid
	gridcpy(grid_new, code_grid);
	pthread_mutex_unlock(&fmutex);
	// END CONCURRENT REGION	
	// reset the new_step_available flag
  if ((error = dlerror()) != NULL) {
    fprintf(stderr, "%s\n", error);
    exit(EXIT_FAILURE);
  }
}

void generator_idle() {
	//pthread_mutex_lock(&gmutex);
	//copy_grid
	gridcpy(grid, grid_new);
	//pthread_mutex_unlock(&gmutex);
	if (gridcmp(grid_new, code_grid) != 0) 
	{	
		cflag++;
		gen_code();
		get_steps(&lin_solve_opt);
	}
}

void start_generator(void *arg) {
	while (1) {  
		generator_idle();
	}
}
