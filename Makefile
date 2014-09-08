all: 
	gcc display.c transformer.c -Wall -lm -lGL -lglut -lGLU -ldl -o fluid
