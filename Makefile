all: 
	gcc display.c transformer.c -Wall -lm -lGL -lglut -lGLU -ldl -fopenmp -O3 -o fluid
