all:
	gcc -std=c99 -o continuum continuum.c -lgsl -lm -lgslcblas
