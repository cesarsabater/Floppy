all: 
	gcc main.c display.c generator.c grid.c simulation_original.c -L/home/cesar/usr/lib -I/home/cesar/usr/include -Wall -lspot -lm -lGL -lglut -lGLU -ldl -losl -lisl -lpthread -fopenmp -O3 -o fluid
	
