#include <stdio.h>
#include <stdlib.h>

extern"C" { double soma_(double *, double *); }

int main(int argc, char * argv[]) {
	double a = atof(argv[1]);
	double b = atof(argv[2]);
	double c = soma_(&a, &b);
	printf("%f", c);
	}