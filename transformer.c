// Compile with -ldl
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dlfcn.h>

#define DENS_FUNCTION "dens_step"
#define VEL_FUNCTION "vel_step"
#define BASE_CODE "simulation_base.c"
#define COMPILED 	"libsimulator_orig.so"

#define max(x,y) (((x) > (y)) ? (x) : (y))

/* extern variables */
extern int N;

extern float ** u, ** v, ** u_prev, ** v_prev;
extern float ** dens, ** dens_prev;
extern float ** vel_max, ** dens_max;

typedef void (*stepfun)(int,float**,float**,float**,float**,float,float);

stepfun get_step(char *fun) {
	void *handle;
  void (*f)(int,float**,float**,float**,float**,float,float);
  char *error;  
	// dynamic compilation of some code (which can be dynamically generated!)
  system("gcc -fPIC -shared -o "COMPILED" "BASE_CODE"");
	// load the dynamic library which have been dynamically generated
  handle = dlopen(COMPILED, RTLD_LAZY);
  if (!handle) {
    fprintf(stderr, "%s\n", dlerror());
    exit(EXIT_FAILURE);
  }
 
  // Get the pointer to the function we want to execute
  f = (void (*)(int,float**,float**,float**,
										float**,float,float))dlsym(handle, fun);
	
  if ((error = dlerror()) != NULL) {
    fprintf(stderr, "%s\n", error);
    exit(EXIT_FAILURE);
  }
  
	return f;
}

stepfun get_dens_step() {
	return get_step(DENS_FUNCTION);
}

stepfun get_vel_step(){
	return get_step(VEL_FUNCTION);
}


