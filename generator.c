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

#define DENS_FUNCTION 	"dens_step"
#define VEL_FUNCTION 		"vel_step"
#define FILENAME     		"solver3.c"
#define BASE_CODE 			"simulation_base.c"
#define ORIGINAL_CODE 	"simulation_original.c"
#define GENERATED_CODE 	"simulation_generated.c"
#define COMPILED 		 		"libsimulator.so"
#define COMPILED_ORIG 	"libsimulator_orig.so"
#define LEVELS 3
#define MAX_ITERATIONS 20
#define FILENAME_LEN 20

/* extern variables */
//main
extern int N;
extern pthread_mutex_t fmutex;
//display
extern float ** u, ** v, ** u_prev, ** v_prev;
extern float ** dens, ** dens_prev;
extern void (*vel_step_opt)(int, float**,float**, float**, float**, float, float);
extern void (*dens_step_opt)(int, float**, float**, float**, float**, float, float);
//grid
extern int G;  // grid size 
extern int slot_size; 
extern int **grid, **grid_aux, **code_grid;	// grid of densities 
extern void refresh_grid(float**);
extern int compare_grids(int **, int **);
extern int iter_from_level(int);


/* local variables */ 
int cflag;
int new_step_available; 
void *handle = NULL;
char obj_name[50];
osl_scop_p scop; 
osl_spot_p spots; 
isl_ctx *ctx;
isl_set **sets; 

// step type
typedef void (*stepfun)(int,float**,float**,float**,float**,float,float);

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
	new_step_available = 0;
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
	for (i = 0; i < LEVELS; i++) {
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
	printf("generacion de codigo!\n");
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

void get_steps(stepfun *d, stepfun *v) {
  char buffer[50]; 
  char *error;
  printf("compilacion!\n");
	// dynamic compilation of some code (which can be dynamically generated!)
	printf("CFLAG: %d\n", cflag);
	
	sprintf(buffer, "sed -i 's/versionXXXX/%d/g' "GENERATED_CODE, cflag);
	system(buffer);
	printf("%s\n", buffer);
	system("gcc -fPIC -shared -o "COMPILED" "GENERATED_CODE);
	printf("gcc -fPIC -shared -o "COMPILED" "GENERATED_CODE"\n");
	
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
	*d = (void (*)(int,float**,float**,float**,
										float**,float,float))dlsym(handle, DENS_FUNCTION);
	*v = (void (*)(int,float**,float**,float**,
										float**,float,float))dlsym(handle, VEL_FUNCTION);
	
	pthread_mutex_unlock(&fmutex);

	// reset the new_step_available flag

  if ((error = dlerror()) != NULL) {
    fprintf(stderr, "%s\n", error);
    exit(EXIT_FAILURE);
  }
}

/*
void get_new_code_if_available() { 
	if (new_code_handle != NULL) {
		printf("charging new code!!!\n");
		//sleep(1);
		if (code_handle != NULL) {
			// references 2 functions
			dlclose(code_handle);
			//dlclose(code_handle);
		} 
		vel_step_opt = vel_step_opt_new;
		dens_step_opt = dens_step_opt_new;
		code_handle = new_code_handle;
		new_code_handle = NULL;
	}
}  
*/

void generator_idle() {
	int i, j;
	if (!compare_grids(grid, code_grid)) { 
		// copy grid
		for (i = 0; i < G; i++)  
		for (j = 0; j < G; j++) {
			// mutex here??
			code_grid[i][j] = grid[i][j];
			// end mutex here??
		}
		//sleep(5);
		gen_code();
		get_steps(&dens_step_opt, &vel_step_opt);
	}
}

void start_generator(void *arg) {
	while (1) {  
		generator_idle();
	}
}