all: 
	gcc display.c simulation_base.c -Wall -lm -lGL -lglut -lGLU -fopenmp -O3 -o fluid
