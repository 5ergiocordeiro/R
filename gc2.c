/*
Solução de sistema de equações lineares (Ax = b) pelo método dos gradientes conjugados.
Emprega Open MPI para aumentar o desempenho.
Uso:
	gcmpi Afile bfile maxiter erro prec
onde
	Afile é o nome do arquivo contendo os valores da matriz A
	bfile é o nome do arquivo contendo os valores do vetor b
	maxiter é o número máximo de iterações permitido
	erro é o erro tolerado
	prec indica se deve ou não usar precondicionamento
		0: não usar precondicionador
		1: usar precondicionador Jacobiano
Limitações:
	a matriz precisa ser simétrica para que o algoritmo funcione,
	para o número de condição foi empregada uma fórmula aproximada,
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mpi.h"


#define bufsize 5000			// para leitura dos dados em arquivo
#define size 37					// tamanho do sistema


int load(char * fname, double * dados, int nrow, int ncol);
double fgetval(char * buffer);


int main(int argc, char * argv[]) {
// Solução de sistema de equações lineares pelo método dos gradientes conjugados
	// Carrega os dados a partir dos arquivos em disco
	double A [size][size], b[size], x[size];				// Ax = b
	double r[size], d[size], q[size], c[size], s[size];		// auxiliares	
	int me, nproc, j, maxiter, prec, mysize;
	double erro;
	MPI_Init(& argc , & argv);
	MPI_Comm_rank(MPI_COMM_WORLD, & me);
	MPI_Comm_size(MPI_COMM_WORLD, & nproc);
	if (me == 0) {
		// ... só o mestre executa a leitura dos dados e dos parâmetros
		int noread;
		noread = load(argv[1], &A[0][0], size, size);
		if (noread > 0) {
			return 1;
			}
		noread = load(argv[2], b, size, 1);
		if (noread > 0) {
			return 1;
			}
		maxiter = atoi(argv[3]) ;
		if (maxiter <= 0) {
			printf("Número máximo de iterações inválido: %d", maxiter);
			return 1;
			}
		erro = atof(argv[4]);
		if (erro <= 0) {
			printf("Erro tolerado inválido: %f", erro);
			return 1;
			}      
		prec = atoi(argv[5]);
		if (prec < 0 || prec > 1) {
			printf("Precondicionador desconhecido: %d", prec);
			return 1;
			}
		// ... faz o tamanho do problema ser divisível pelo número de processadores
		my_size = size + size % (nproc - 1);
		// ... comunica os valores aos demais
		MPI_Bcast (& my_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (& prec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (& erro, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast (& maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (& A [0][0] , my_size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast (b , my_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
	else {
		// ... os demais processadores recebem os valores passados pelo mestre	
		MPI_Recv (& my_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Recv (& prec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Recv (& erro, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Recv (& maxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Recv (& A [0][0] , my_size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Recv (b , my_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
	// Resolve o sistema
	// ... inicialização
	int i, k;
	double sum, p, pt;
	for (i = me * my_size; i < (me + 1) * my_size ; ++ i) {
		// ... calcula o precondicionador
		if (prec == 1) {
			// precondicionador Jacobiano
			c[i] = 1.0 / A[i][i] ;
			}
		// ... inicializa o vetor solução (x0)
		// (teoricamente, pode ser qualquer coisa)
		x[i] = 0;
		}
	// ... todos anunciam os vetores inicializados
	MPI_Bcast (c + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);
	MPI_Bcast (x + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);
	// ... todos recebem os vetores inicializados
	MPI_Reduce(c, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	MPI_Reduce(x, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	pt = p = 0;
	for (i = me * my_size; i < (me + 1) * my_size ; ++ i) {
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
			pt += r[i] * r[i];
			}
		if (prec == 1) {
			// d(0) = c' r(0)
			// p(0) = r(0)' d(0)
			d[i] = c[i] * r[i];
			pt += r[i] * d[i];
			}
		}
	// ... todos anunciam os valores calculados
	MPI_Bcast(r + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);		
	MPI_Bcast(d + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);
	MPI_Bcast(& pt, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);		
	// ... todos recebem os valores calculados
	MPI_Reduce(r, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	MPI_Reduce(d, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	MPI_Reduce(& pt, & p, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
	// ... executa a iteração
	for (j = 1; j < maxiter && fabs(p) > erro; ++ j) {
		// q(j) = A d(j)
		// qq(j) = d(j)' A d(j)
		double qq = qt = 0 ; 
		for (i = me * my_size; i < (me + 1) * my_size ; ++ i) {
			sum = 0;
			for (k = 0; k < size; ++ k) {
				sum = sum + A[i][k] * d[k];
				}
			qq += d[i] * sum;
			q[i] = sum;
			}
		// ... todos anunciam os valores calculados
		MPI_Bcast(q + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);		
		MPI_Bcast(& qq, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);		
		// ... todos recebem os valores calculados
		MPI_Reduce(q, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		MPI_Reduce(& qt, & qq, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );	
		// x(j+1) = x(j) + p/qq d(j)
		// r(j+1) = r(i) - p/qq A d(j)
		double alfa = p/qt;
		for (i = me * my_size; i < (me + 1) * my_size ; ++ i) {
			x[i] += alfa * d[i];
			r[i] -= alfa * q[i];
			}
		// ... todos anunciam os valores calculados
		MPI_Bcast(x + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);		
		MPI_Bcast(r + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);		
		// ... todos recebem os valores calculados
		MPI_Reduce(x, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		MPI_Reduce(r, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		double lastp = p;
		pt = p = 0;
		for (i = me * my_size; i < (me + 1) * my_size ; ++ i) {
			if (prec == 0) {
				// p(j+1) = r(j+1)' r(j+1)
				pt += r[i] * r[i];
				}
			if (prec == 1) {
				// s(j+1) = c' r(j+1)
				// p(j+1) = r(j+1)' s(j+1)
				s[i] = c[i] * r[i];
				pt += r[i] * s[i];
				}
			}
		// ... todos anunciam os valores calculados
		if (prec == 1) {
			MPI_Bcast(s + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);
			}
		MPI_Bcast(& pt, 1, MPI_DOUBLE, me, MPI_COMM_WORLD);		
		// ... todos recebem os valores calculados
		if (prec == 1) {
			MPI_Reduce(s, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			}
		MPI_Reduce(& pt, & p, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );	
		double beta = p/lastp;
		for (i = me * my_size; i < (me + 1) * my_size ; ++ i) {
			if (prec == 0) {
				// d(j+1) = r(j+1) + p(j+1)/p(j) d(j)
				d[i] = r[i] + beta * d[i];
				}
			if (prec == 1) {
				// d(j+1) = s(j+1) + p(j+1)/p(j) d(j)
				d[i] = s[i] + beta * d[i];
				}
			}
		// ... todos anunciam os valores calculados
		MPI_Bcast(d + me * my_size, my_size, MPI_DOUBLE, me, MPI_COMM_WORLD);		
		// ... todos recebem os valores calculados
		MPI_Reduce(d, 0, my_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		}
	if (me != 0) {
		return 0;
		}
	// ... o mestre exibe a solução
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
	printf("] \nErro < %f, ka = %f", fabs(p), ka);
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
