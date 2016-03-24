/*
Testa a implementação de funções de diversos níveis da BLAS.
Uso:
	myblas rep n size1 size2 ... sizen
onde
	rep é o número de repetições a executar
	n é o nível BLAS, um inteiro de 1 a 3
	size1, size2, sizen são o número de elementos nos vetores e matrizes de teste
O programa calcula o número de operações de ponto flutuante teoricamente executadas.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <windows.h>


// #define IMPRIMIR


void Init (int size, double * vals) {
// Inicializa um conjunto de valores aleatoriamente
	while (size -- > 0)
		* vals ++ = rand() % 10;
	}

double ddot(double * x, double * y, int size) {
// Retorna o produto escalar de dois vetores
	double sum = 0;
	while (size -- > 0)
		sum += (* x ++) * (* y ++);
	return sum;
	}
	
void dgemv (double * A, double * x, double * y, int nrows, int ncols) {
// Calcula o produto de uma matriz por um vetor
	while (nrows -- > 0) {
		* y ++ = ddot(A, x, ncols);
		A += ncols;
		}
	return ;
	}
	
void dgemm (double * A, double * B, double * C, int nrowsA, int nrowsB, int ncolsC) {
// Calcula o produto de duas matrizes
	double * ptB, * tB = (double *) malloc(nrowsB * ncolsC * sizeof(double));
	if (tB == NULL) {
		printf("Falha ao alocar a memória de trabalho necessária");
		exit(1);			
		}
	ptB = tB;
	for (int i = 0 ; i < nrowsB ; ++ i) {
		for (int j = 0 ; j < ncolsC ; ++ j) {
			* ptB ++ = * (B + j * ncolsC + i);
			}
		}
	while (nrowsA -- > 0) {
		ptB = tB;
		for (int i = 0; i < nrowsB; ++ i) {
			* C ++ = ddot(A, ptB, ncolsC);
			ptB += ncolsC;
			}
		A += nrowsB;
		}
	}
	
int main(int argc, char * argv[]) {
// Testa a implementação de funções de diversos níveis da BLAS
	int level = atoi(argv[1]);
	if (level < 1 || level > 3) {
		printf("Nível inválido: %d", level);
		return 1;
		}
	int rep = atoi(argv[2]);
	if (rep <= 0) {
		printf("Número de repetições inválido: %d", rep);
		return 1;
		}
    HANDLE hProcess = GetCurrentProcess();
    FILETIME ftCreation, ftExit, ftKernel, ftUser1, ftUser2;
    SYSTEMTIME stUser1,stUser2;
	double gflop;
	if (level == 1) {
		if (argc != 4 ) {
			printf("Número inválido de argumentos");
			return 1;
			}
		int ncols = atoi(argv[3]);
		gflop = rep * ncols/1e9;
		int mcols = ncols/10;
		double * x = (double *) malloc(ncols * sizeof(double));
		double * y = (double *) malloc(ncols * sizeof(double));
		if (x == NULL || y == NULL) {
			printf("Falha ao alocar a memória de trabalho necessária");
			return 1;			
			}
		Init(ncols, x);
		Init(ncols, y);
		GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser1);
		for (int i = 0; i < rep; ++ i)
			double result = ddot(x, y, ncols);
		GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser2);
#ifdef IMPRIMIR
		printf("[");
		for (int i = 0; i < ncols; ++ i) {
			if (i > 0) printf(" ");
			printf("%f" , x[i]);
			}
		printf("] .\n[");
		for (int i = 0; i < ncols; ++ i) {
			if (i > 0) printf(" ");
			printf("%f" , y[i]);
			}
		printf("]\n\t= %f", result);
#endif
		}
	if (level == 2) {
		if (argc != 5) {
			printf("Número inválido de argumentos");
			return 1;
			}
		int nrows = atoi(argv[3]);
		int ncols = atoi(argv[4]);
		gflop = rep * (2.0 * ncols - 1) * nrows / 1e9;
		double * pA, * A = (double *) malloc(nrows * ncols * sizeof(double));
		double * px, * x = (double *) malloc(ncols * sizeof(double));
		double * py, * y = (double *) malloc(nrows * sizeof(double));
		if (A == NULL || x == NULL || y == NULL) {
			printf("Falha ao alocar a memória de trabalho necessária");
			return 1;			
			}
		Init(nrows * ncols, A);
		Init(ncols, x);
		GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser1);
		for (int i = 0; i < rep; ++ i)
			dgemv(A, x, y, nrows, ncols);
		GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser2);
#ifdef IMPRIMIR
		pA = A;
		for (int i = 0; i < nrows; ++ i) {
			printf("|");
			for (int j = 0; j < ncols; ++ j) {
				if (j > 0) printf(" ");
				printf("%f" , * pA ++) ;
				}
			printf("|\n");
			}
		printf("x [");
		px = x;
		for (int i = 0; i < ncols; ++ i) {
			if (i > 0) printf(" ");
			printf("%f" , * px ++);
			}
		printf("]\n\t= [");
		py = y;
		for (int i = 0; i < nrows; ++ i) {
			if (i > 0) printf(" ");
			printf("%f" , * py ++);
			}
		printf("]\n");
#endif
		}
	if (level == 3) {
		if (argc != 6) {
			printf("Número inválido de argumentos");
			return 1;
			}
		int nrowsA = atoi(argv[3]);
		int nrowsB = atoi(argv[4]);
		int ncolsC = atoi(argv[5]);
		gflop = rep * (2.0 * nrowsB - 1) * nrowsA * ncolsC / 1e9;
		double * pA, * A = (double *) malloc(nrowsA * nrowsB * sizeof(double));
		double * pB, * B = (double *) malloc(nrowsB * ncolsC * sizeof(double));
		double * pC, * C = (double *) malloc(nrowsA * ncolsC * sizeof(double));
		if (A == NULL || B == NULL || C == NULL) {
			printf("Falha ao alocar a memória de trabalho necessária");
			return 1;			
			}
		Init(nrowsA * nrowsB, A);
		Init(nrowsB * ncolsC, B);
		GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser1);
		for (int i = 0; i < rep; ++ i)
			dgemm(A, B, C, nrowsA, nrowsB, ncolsC);
		GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser2);
#ifdef IMPRIMIR
		pA = A;
		for (int i = 0; i < nrowsA; ++ i) {
			printf("|");
			for (int j = 0; j < nrowsB; ++ j) {
				if (j > 0) printf(" ");
				printf("%f" , * pA ++);
				}
			printf("|\n");
			}
		printf("x\n");
		pB = B;
		for (int i = 0; i < nrowsB; ++ i) {
			printf("|");
			for (int j = 0; j < ncolsC; ++ j) {
				if (j > 0) printf(" ");
				printf("%f" , * pB ++);
				}
			printf("|\n");
			}
		printf("=\n");
		pC = C;
		for (int i = 0; i < nrowsA; ++ i) {
			printf("|");
			for (int j = 0; j < ncolsC; ++ j) {
				if (j > 0) printf(" ");
				printf("%f" , * pC ++);
				}
			printf("|\n");
			}
#endif
		}
	FileTimeToSystemTime(& ftUser1, & stUser1);
	FileTimeToSystemTime(& ftUser2, & stUser2);
	printf("GFlops = %f, User time = %d ms \n", gflop,
		stUser2.wSecond * 1000 - stUser1.wSecond * 1000 + stUser2.wMilliseconds - stUser1.wMilliseconds);
	return 0;
	}
