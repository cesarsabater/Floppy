all: 
	gcc display.c transformer.c -L/home/cesar/usr/lib -I/home/cesar/usr/include -Wall -lspot -lm -lGL -lglut -lGLU -ldl -losl -lisl -o fluid
	

