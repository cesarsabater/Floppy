# Compiller for Lin Solv
RUNTIMECOMPILER=gcc
RUNTIMEOPTS=-O3
# ACR Compiller
ACRCOMPILER=gcc
ACROPTS=-O3
# Trivial Static Optimization Compiller
STATICOPTCC=gcc
SOPTIMOPT=-O3

all:
	$(RUNTIMECOMPILER) $(RUNTIMEOPTS) -c lin_solve_generated.c
	$(ACRCOMPILER) $(ACROPTS) controller.c simulation_opt.c -fopenmp -o sim_opt lin_solve_generated.o
	$(STATICOPTCC) $(OPTIMOPT) controller.c simulation_trivial.c -fopenmp -o sim_triv

## old stuff
#all:
#	tcc generator.c -L/home/cesar/usr/lib -I/home/cesar/usr/include -Wall -lspot -lm -losl -lisl -fopenmp -o gen
