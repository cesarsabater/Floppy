#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <isl/ctx.h>
#include <isl/options.h>
#include <isl/printer.h>
#include <isl/set.h>

#include <spot/spot.h>
#include <osl/osl.h>
#include <osl/extensions/doi.h>

#define BASE_CODE 			"lin_solve_base.c"
#define GENERATED_CODE 	"lin_solve_generated.c"

/* global variables */

//grid
int N = 200;
int G = 10;  // grid size 
int slot_size; 
int **grid;
//polyhedral data structures
osl_scop_p scop; 
osl_spot_p spots; 
isl_ctx *ctx;
isl_set **sets; 

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
}

void grid_init()
{
	slot_size	= N / G;
	alloc_data();
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

int get_levels() 
{
	return 3;
}

int max_outer_loop_iterations() {
	return 20;
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
	// TODO: FREE SETS!!!
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

int main() { 
	grid_init();
	read_grid();
	gen_code();
	return 0;
}
