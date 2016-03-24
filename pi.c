/*
Cálculo do valor de PI por integração trapezoidal de função de referência
Uso:
	pi n
onde n é o número de divisões do eixo x. O valor default é 400, que dá o valor exato.
Pode usar OpenMP para aumentar o desempenho.
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main(int argc, char * argv[]) {
// Calcula o valor de PI por integração trapezoidal de função de referência
	int divisions;
	if (argc == 1) {
		divisions = 400;
		}
	else {
		divisions = atol(argv[1]);
		}
	long i;
	double sum = 0, h = 1.0 / divisions;
	// Para desativar a paralelização por OpenMp, remover a linha abaixo
	#pragma omp parallel for private(i) reduction(+:sum)	
	for(i = 1; i <= divisions; ++i) {
		sum = sum + 1 / (1 + i * h * i * h) + 1 / (1 + (i-1) * h * (i-1) * h );
		}
	sum *= 2.0/divisions;
	printf("PI = %f \n",sum);
	}	