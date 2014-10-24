#options
CC = gcc
OPT = -O3

all: 
	$(CC) display.c simulation_generated.c grid.c -Wall -lm -lGL -lglut -lGLU -fopenmp $(OPT) -o fluid
