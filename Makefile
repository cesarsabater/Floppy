all: 
	gcc display.c simulation_generated.c grid.c -Wall -lm -lGL -lglut -lGLU -fopenmp -O3 -o fluid
