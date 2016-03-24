/*
Solução de sistema de equações lineares (Ax = b) pelo método dos gradientes conjugados.
Pode empregar OpenMP para aumentar o desempenho.
Uso:
	gc Afile bfile maxiter erro prec omp
onde
	Afile é o nome do arquivo contendo os valores da matriz A
	bfile é o nome do arquivo contendo os valores do vetor b
	maxiter é o número máximo de iterações permitido
	erro é o erro tolerado
	prec indica se deve ou não usar precondicionamento
		0: não usar precondicionador
		1: usar precondicionador Jacobiano
	omp indica se vai usar o OpenMP
		0: não usar
		1: usar
Limitações:
	a matriz precisa ser simétrica para que o algoritmo funcione,
	para o número de condição foi empregada uma fórmula aproximada,
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>

#include <sys/time.h>
#include <windows.h>
#include <time.h>

#define bufsize 5000			// para leitura dos dados em arquivo
#define size 37					// tamanho do sistema


int load(char * fname, double * dados, int nrow, int ncol);
double fgetval(char * buffer);


int main(int argc, char * argv[]) {
// Solução de sistema de equações lineares pelo método dos gradientes conjugados
	// Carrega os dados a partir dos arquivos em disco
	double wstart, wend;
	
	double A [size][size], b[size], x[size];				// Ax = b
	double r[size], d[size], q[size], c[size], s[size];		// auxiliares	
	int noread;
	noread = load(argv[1], &A[0][0], size, size);
	if (noread > 0) {
		return 1;
		}
	noread = load(argv[2], b, size, 1);
	if (noread > 0) {
		return 1;
		}
	int j, maxiter = atoi(argv[3]) ;
	if (maxiter <= 0) {
		printf("Número máximo de iterações inválido: %d", maxiter);
		return 1;
		}
	double erro = atof(argv[4]);
	if (erro <= 0) {
		printf("Erro tolerado inválido: %f", erro);
		return 1;
		}
	int prec = atoi(argv[5]);
	if (prec < 0 || prec > 1) {
		printf("Precondicionador desconhecido: %d", prec);
		return 1;
		}
	int omp = atoi(argv[6]);
	if (omp == 0) {
		omp_set_num_threads(1);
		}	
	// Resolve o sistema
    HANDLE hProcess = GetCurrentProcess();
    FILETIME ftCreation, ftExit, ftKernel, ftUser1, ftUser2;
    SYSTEMTIME stUser1,stUser2;
    GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser1);
	wstart = omp_get_wtime();
	// ... inicialização
	int i, k, jjj = 0;
	double sum, p, qq;
	for(int count = 0; count < 1; ++ count) {
	#pragma omp parallel for private(i)
	for (i = 0; i < size ; ++ i) {
		// ... calcula o precondicionador
		if (prec == 1) {
			// precondicionador Jacobiano
			c[i] = 1.0 / A[i][i] ;
			}
		// ... inicializa o vetor solução (x0)
		// (teoricamente, pode ser qualquer coisa)
		x[i] = 0;
		}
	p = 0;
	#pragma omp parallel for private(i,sum,k) reduction(+:p)
	for (i = 0; i < size ; ++ i) {
		// ... calcula os vetores iniciais r(0) e d(0)
		sum = 0;
		// r(0) = b - A x(0)
		for (k = 0; k < size; ++ k) {
			sum = sum + A[i][k] * x[k];
			}
		r[i] = b[i] - sum;
		if (prec == 0) {
			// d(0) = r(0)
			// p(0) = r(0)' r(0)
			d[i] = r[i];
			p += r[i] * r[i];
			}
		if (prec == 1) {
			// d(0) = c' r(0)
			// p(0) = r(0)' d(0)
			d[i] = c[i] * r[i];
			p += r[i] * d[i];
			}
		}
	// ... executa a iteração
	for (j = 1; j < maxiter && fabs(p) > erro; ++ j) {
		// q(j) = A d(j)
		// qq(j) = d(j)' A d(j)
		qq = 0; 
		#pragma omp parallel for private(i,sum,k) reduction(+:qq)
		for (i = 0; i < size; ++ i) {
			sum = 0;
			for (k = 0; k < size; ++ k) {
				sum = sum + A[i][k] * d[k];
				}
			qq += d[i] * sum;
			q[i] = sum;
			}
		// x(j+1) = x(j) + p/qq d(j)
		// r(j+1) = r(i) - p/qq A d(j)
		double alfa = p/qq;
		#pragma omp parallel for private(i)
		for (i = 0; i < size; ++ i) {
			x[i] += alfa * d[i];
			r[i] -= alfa * q[i];
			}
		double lastp = p;
		p = 0;
		#pragma omp parallel for private(i) reduction(+:p)
		for (i = 0; i < size; ++ i) {
			if (prec == 0) {
				// p(j+1) = r(j+1)' r(j+1)
				p += r[i] * r[i];
				}
			if (prec == 1) {
				// s(j+1) = c' r(j+1)
				// p(j+1) = r(j+1)' s(j+1)
				s[i] = c[i] * r[i];
				p += r[i] * s[i];
				}
			}
		double beta = p/lastp;
		#pragma omp parallel for private(i)
		for (i = 0; i < size; ++ i) {
			if (prec == 0) {
				// d(j+1) = r(j+1) + p(j+1)/p(j) d(j)
				d[i] = r[i] + beta * d[i];
				}
			if (prec == 1) {
				// d(j+1) = s(j+1) + p(j+1)/p(j) d(j)
				d[i] = s[i] + beta * d[i];
				}
			}
		// if (omp_get_thread_num() == 0) printf("%Iteração %d %d: erro = %f \n", count, j, p);
		}
	}
    GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser2);	
	wend = omp_get_wtime();
    FileTimeToSystemTime(& ftUser1, & stUser1);
    FileTimeToSystemTime(& ftUser2, & stUser2);
	double twall = 1000.0 * (wend - wstart);
    double tuser = 1000.0 * (stUser2.wSecond - stUser1.wSecond) + stUser2.wMilliseconds - stUser1.wMilliseconds;
	printf("tempo gasto = %f ms, user time = %f ms %d\n", twall, tuser, jjj);
	// Exibe a solução
	double ka, val, max = 0, min = 1e9;
	printf("Solução após %d iterações: [", j - 1);
	for (i = 0; i < size ; ++ i) {
		if (i > 0) {
			printf(", ");
			}
		printf("%f", x[i]);
		if (prec == 0) {
			val = fabs(A[i][i]);
			if (val > max) max = val;
			if (val < min) min = val;
			}
		}
	if (prec == 0) {
		ka = max / min;
		}
	if (prec == 1) {
		ka = 1;
		}
	printf("] \nErro < %f, ka = %f, ", fabs(p), ka);
	return 0;
	}
	
int load(char * fname, double * dados, int nrow, int ncol) {
// Carrega os dados do arquivo "fname" na matriz "dados", de dimensão "nrow" x "ncol"
	char * pchar, buf [bufsize];
	FILE * fp = fopen (fname, "r");
	if (fp == NULL) {
		printf("Não conseguiu ler o arquivo %s. \n", fname);
		return 1;
		}
	fgets(buf,bufsize,fp);
	fgets(buf,bufsize,fp);
	for (int i = 0; i < nrow; ++ i) {
		fgets(buf,bufsize,fp);
		for (int k = 0; k < ncol; ++ k) {
			* dados = fgetval(buf);
			dados ++ ;
			}
		}
	fclose(fp);
	return 0;
	}


double fgetval(char * buffer) {
// Extrai o primeiro valor existente no "buffer"
	char * pchar = buffer, valor [bufsize], * pvalor = valor;
	while (isblank(*pchar) ) {
		pchar ++;
		}
	while (! isblank(*pchar) )	{
		* pvalor ++ = * pchar ++;
		}
	* pvalor = '\0';
	strcpy(buffer,pchar);
	return atof(valor);
	}
