/*
solvegauss.c
Lê um sistema gerado pelo MATLAB, resolve-o pelo método de eliminação de Gauss e calcula o determiante e a norma 2 do resultado em diversas precisões.
Uso:
	solvegauss
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <quadmath.h>

#define bufsize 10000			// para leitura dos dados em arquivo
#define POSROWNBR 7				// posição do número de linhas no arquivo
#define POSCOLNBR 10			// posição do número de colunas no arquivo


double fgetval(char ** pbuffer);
double dmnorm2(double * pmat, int nrow, int ncol);
float fmnorm2(float * pmat, int nrow, int ncol);
long double ldmnorm2(long double * pmat, int nrow, int ncol);
float * fmcopy(double * psrc, int nrows, int ncols);
float * fmtrisolve(float * pmat, int nrows, int ncols, bool superior);
double * dmtrisolve(double * pmat, int nrows, int ncols, bool superior);
long double * ldmtrisolve(long double * pmat, int nrows, int ncols, bool superior);
long double * ldmcopy(double * psrc, int nrows, int ncols);
double * lermat(const char * fname, int * nrows, int * ncolA);
int main(int argc, const char * argv[]);
float * fsolveG(float * psrc, int rank, float * pdet);
double * dsolveG(double * psrc, int rank, double * pdet);
long double * ldsolveG(long double * psrc, int rank, long double * pdet);
int ffindmax(float * pmat, int nrows, int ncols, int pos, bool colmode, int start);	
int dfindmax(double * pmat, int nrows, int ncols, int pos, bool colmode, int start);	
int ldfindmax(long double * pmat, int nrows, int ncols, int pos, bool colmode, int start);	
void fchangerows(float * pmat, int rows, int ncols, int row1, int row2);
void dchangerows(double * pmat, int rows, int ncols, int row1, int row2);
void ldchangerows(long double * pmat, int rows, int ncols, int row1, int row2);
float * f2tri(float * psrc, int rank, float * pdet);


static int debuglevel_ = 0;

int main(int argc, const char * argv[]) {
	// Lê o sistema de entrada
	float * pAf, * pCf, fnorm;
	double * pAd, * pCd, dnorm;
	long double * pAld, * pCld, ldnorm;
	// __float128 * pAq, * pBq, * pCq, qnorm;
	int nrowA, ncolA;
	pAd = lermat("S", & nrowA, & ncolA);
	if (pAd == NULL) {
		return 1;
		}
	// Verifica se pode ser resolvido
	if (ncolA != nrowA + 1) {
		printf("O sistema não podem ser resolvido, porque as dimensões são incompatíveis: (%d,%d) \n", nrowA, ncolA);
		return 2;
		}
	// Cria versões em diversas precisões
	pAf = fmcopy(pAd, nrowA, ncolA);
	pAld = ldmcopy(pAd, nrowA, ncolA);
	// Resolve o sistema
	float fdet;
	double ddet;
	long double lddet;	
	pCf = fsolveG(pAf, nrowA, & fdet);
	pCd = dsolveG(pAd, nrowA, & ddet);
	pCld = ldsolveG(pAld, nrowA, & lddet);
	// Calcula a norma 2 dos resultados
	fnorm = fmnorm2(pCf, nrowA, 1);
	dnorm = dmnorm2(pCd, nrowA, 1);
	ldnorm = ldmnorm2(pCld, nrowA, 1);
	printf("Determinante da matriz: 32 bits = %f, 64 bits = %f, 80 bits = %Lf \n", fdet, ddet, lddet);	
	printf("Norma 2 do resultado: 32 bits = %f, 64 bits = %f, 80 bits = %Lf \n", fnorm, dnorm, ldnorm);	
	return 0;
	}

// Funções para cópia das matrizes em diversas precisões	
float * fmcopy(double * psrc, int nrows, int ncols) {
// Retorna uma cópia em precisão simples (32 bits) da matriz 'psrc'
	float * pdst, * result;
	int size = nrows * ncols;
	result = pdst = (float *) malloc (size * sizeof(float));
	if (pdst == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", nrows, ncols);
		return NULL;
		}
	while (size -- > 0) {
		* pdst ++ = * psrc ++;
		}
	return result;
	}
	
long double * ldmcopy(double * psrc, int nrows, int ncols) {
// Retorna uma cópia em precisão estendida (80 bits) da matriz 'psrc'
	long double * pdst, * result;
	int size = nrows * ncols;
	result = pdst = (long double *) malloc (size * sizeof(long double));
	if (pdst == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", nrows, ncols);
		return NULL;
		}
	while (size -- > 0) {
		double dvalor;
		long double ldvalor;
		dvalor = * psrc ++;
		ldvalor = dvalor;
		* pdst ++ = ldvalor;
		}
	return result;
	}

// Funções para solução dos sistemas em diversas precisões
float * f2tri(float * psrc, int rank, float * pdet) {
	int ncols = rank + 1;
	float * pval = (float *) calloc(rank * ncols, sizeof(float));
	if (pval == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", rank, ncols);
		return NULL;
		}
	for (int i = 0; i < rank; ++ i) {
		for (int j = 0; j <= rank; ++ j) {
			pval[i * ncols + j] = psrc[i * (rank + 1) + j];
			}
		}
	if (debuglevel_ >= 2) {
		printf("Inicialização \n");
		for (int i = 0; i < rank; ++ i) {
			for (int j = 0; j < ncols; ++ j) {
				printf(" %f ", pval[i * ncols + j]);
				}
			printf("\n");
			}
		}
	bool sinal = false;
	float det = 1, maxval;
	for (int j = 0; j < rank - 1; ++ j) {
		int maxrow = ffindmax(pval, rank, ncols, j, true, j);
		maxval = pval[maxrow * ncols + j];
		if (maxval == 0) {
			printf("Matriz singular! \n");
			return NULL;
			}
		if (j != maxrow) {
			sinal = ! sinal;
			fchangerows(pval, rank, ncols, j, maxrow);
			if (debuglevel_ >= 2) {
				printf("Coluna %d pivoteamento\n", j);
				for (int i = 0; i < rank; ++ i) {
					for (int k = 0; k < ncols; ++ k) {
						printf(" %f ", pval[i * ncols + k]);
						}
					printf("\n");
					}
				}
			}
		for (int i = j + 1; i < rank; ++ i) {
			double multiplier = pval[i * ncols + j] / maxval;
			pval[i * ncols + j] = 0;
			for (int k = j + 1; k <= rank; ++ k) {
				pval[i * ncols + k] -= pval[j * ncols + k] * multiplier;
				}
			}
		if (debuglevel_ >= 2) {
			printf("Coluna %d eliminação. Pivô = %f\n", j, maxval);
			for (int i = 0; i < rank; ++ i) {
				for (int k = 0; k < ncols; ++ k) {
					printf(" %f ", pval[i * ncols + k]);
					}
				printf("\n");
				}
			}
		}
	for (int i = 0; i < rank; ++ i) {
		det *= pval[i * ncols + i];
		}
	* pdet = det * (sinal ? -1 : 1);
	return pval;
	}

float * fsolveG(float * psrc, int rank, float * pdet) {
	float * pTS = f2tri(psrc, rank, pdet);
	if (pTS == NULL) {
		return NULL;
		}
	float * result = (float *) calloc(rank, sizeof(float));
	if (result == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x 1! \n", rank);
		return NULL;
		}
	result = fmtrisolve(pTS, rank, rank + 1, true);
	if (debuglevel_ >= 2) {
		printf("Resultado \n");
		for (int i = 0; i < rank; ++ i) {
			printf(" %f ", result[i]);
			}
		printf("\n");
		}
	return result;
	}

double * d2tri(double * psrc, int rank, double * pdet) {
	int ncols = rank + 1;
	double * pval = (double *) calloc(rank * ncols, sizeof(double));
	if (pval == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", rank, ncols);
		return NULL;
		}
	for (int i = 0; i < rank; ++ i) {
		for (int j = 0; j <= rank; ++ j) {
			pval[i * ncols + j] = psrc[i * (rank + 1) + j];
			}
		}
	if (debuglevel_ >= 2) {
		printf("Inicialização \n");
		for (int i = 0; i < rank; ++ i) {
			for (int j = 0; j < ncols; ++ j) {
				printf(" %f ", pval[i * ncols + j]);
				}
			printf("\n");
			}
		}
	bool sinal = false;
	double det = 1, maxval;
	for (int j = 0; j < rank - 1; ++ j) {
		int maxrow = dfindmax(pval, rank, ncols, j, true, j);
		maxval = pval[maxrow * ncols + j];
		if (maxval == 0) {
			printf("Matriz singular! \n");
			return NULL;
			}
		if (j != maxrow) {
			sinal = ! sinal;
			dchangerows(pval, rank, ncols, j, maxrow);
			if (debuglevel_ >= 2) {
				printf("Coluna %d pivoteamento\n", j);
				for (int i = 0; i < rank; ++ i) {
					for (int k = 0; k < ncols; ++ k) {
						printf(" %f ", pval[i * ncols + k]);
						}
					printf("\n");
					}
				}
			}
		for (int i = j + 1; i < rank; ++ i) {
			double multiplier = pval[i * ncols + j] / maxval;
			pval[i * ncols + j] = 0;
			for (int k = j + 1; k <= rank; ++ k) {
				pval[i * ncols + k] -= pval[j * ncols + k] * multiplier;
				}
			}
		if (debuglevel_ >= 2) {
			printf("Coluna %d eliminação. Pivô = %f\n", j, maxval);
			for (int i = 0; i < rank; ++ i) {
				for (int k = 0; k < ncols; ++ k) {
					printf(" %f ", pval[i * ncols + k]);
					}
				printf("\n");
				}
			}
		}
	for (int i = 0; i < rank; ++ i) {
		det *= pval[i * ncols + i];
		}
	* pdet = det * (sinal ? -1 : 1);
	return pval;
	}

double * dsolveG(double * psrc, int rank, double * pdet) {
	double * pTS = d2tri(psrc, rank, pdet);
	if (pTS == NULL) {
		return NULL;
		}
	double * result = (double *) calloc(rank, sizeof(double));
	if (result == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x 1! \n", rank);
		return NULL;
		}
	result = dmtrisolve(pTS, rank, rank + 1, true);
	if (debuglevel_ >= 2) {
		printf("Resultado \n");
		for (int i = 0; i < rank; ++ i) {
			printf(" % f ", result[i]);
			}
		printf("\n");
		}
	return result;
	}

long double * ld2tri(long double * psrc, int rank, long double * pdet) {
	int ncols = rank + 1;
	long double * pval = (long double *) calloc(rank * ncols, sizeof(long double));
	if (pval == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", rank, ncols);
		return NULL;
		}
	for (int i = 0; i < rank; ++ i) {
		for (int j = 0; j <= rank; ++ j) {
			pval[i * ncols + j] = psrc[i * (rank + 1) + j];
			}
		}
	if (debuglevel_ >= 2) {
		printf("Inicialização \n");
		for (int i = 0; i < rank; ++ i) {
			for (int j = 0; j < ncols; ++ j) {
				printf(" %Lf ", pval[i * ncols + j]);
				}
			printf("\n");
			}
		}
	bool sinal = false;
	long double det = 1, maxval;
	for (int j = 0; j < rank - 1; ++ j) {
		int maxrow = ldfindmax(pval, rank, ncols, j, true, j);
		maxval = pval[maxrow * ncols + j];
		if (maxval == 0) {
			printf("Matriz singular! \n");
			return NULL;
			}
		if (j != maxrow) {
			sinal = ! sinal;
			ldchangerows(pval, rank, ncols, j, maxrow);
			if (debuglevel_ >= 2) {
				printf("Coluna %d pivoteamento\n", j);
				for (int i = 0; i < rank; ++ i) {
					for (int k = 0; k < ncols; ++ k) {
						printf(" %Lf ", pval[i * ncols + k]);
						}
					printf("\n");
					}
				}
			}
		for (int i = j + 1; i < rank; ++ i) {
			long double multiplier = pval[i * ncols + j] / maxval;
			pval[i * ncols + j] = 0;
			for (int k = j + 1; k <= rank; ++ k) {
				pval[i * ncols + k] -= pval[j * ncols + k] * multiplier;
				}
			}
		if (debuglevel_ >= 2) {
			printf("Coluna %d eliminação. Pivô = %Lf\n", j, maxval);
			for (int i = 0; i < rank; ++ i) {
				for (int k = 0; k < ncols; ++ k) {
					printf(" %Lf ", pval[i * ncols + k]);
					}
				printf("\n");
				}
			}
		}
	for (int i = 0; i < rank; ++ i) {
		det *= pval[i * ncols + i];
		}
	* pdet = det * (sinal ? -1 : 1);
	return pval;
	}

long double * ldsolveG(long double * psrc, int rank, long double * pdet) {
	long double * pTS = ld2tri(psrc, rank, pdet);
	if (pTS == NULL) {
		return NULL;
		}
	long double * result = (long double *) calloc(rank, sizeof(long double));
	if (result == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x 1! \n", rank);
		return NULL;
		}
	result = ldmtrisolve(pTS, rank, rank + 1, true);
	if (debuglevel_ >= 2) {
		printf("Resultado \n");
		for (int i = 0; i < rank; ++ i) {
			printf(" %Lf ", result[i]);
			}
		printf("\n");
		}
	return result;
	}

int ffindmax(float * pmat, int nrows, int ncols, int pos, bool colmode, int start) {
	float maxval = 0;
	int i, address, maxpos = start;
	int size = colmode ? nrows : ncols;
	for (i = start; i < size; ++ i) {
		if (colmode) {
			address = i * ncols + pos;
			}
		else {
			address = pos * nrows + i;
			}
		float value = fabs(pmat[address]);
		if (maxval < value) {
			maxval = value;
			maxpos = i;
			}
		}
	return maxpos;
	}
	
void fchangerows(float * pmat, int rows, int ncols, int row1, int row2) {
	int row1pos = row1 * ncols;
	int row2pos = row2 * ncols;
	for (int j = 0; j < ncols; ++ j) {
		float value = pmat[row1pos + j];
		pmat[row1pos + j] = pmat[row2pos + j];
		pmat[row2pos + j] = value;
		}
	}

int dfindmax(double * pmat, int nrows, int ncols, int pos, bool colmode, int start) {
	double maxval = 0;
	int i, address, maxpos = start;
	int size = colmode ? nrows : ncols;
	for (i = start; i < size; ++ i) {
		if (colmode) {
			address = i * ncols + pos;
			}
		else {
			address = pos * nrows + i;
			}
		double value = fabs(pmat[address]);
		if (maxval < value) {
			maxval = value;
			maxpos = i;
			}
		}
	return maxpos;
	}
	
void dchangerows(double * pmat, int rows, int ncols, int row1, int row2) {
	int row1pos = row1 * ncols;
	int row2pos = row2 * ncols;
	for (int j = 0; j < ncols; ++ j) {
		double value = pmat[row1pos + j];
		pmat[row1pos + j] = pmat[row2pos + j];
		pmat[row2pos + j] = value;
		}
	}

int ldfindmax(long double * pmat, int nrows, int ncols, int pos, bool colmode, int start) {
	long double maxval = 0;
	int i, address, maxpos = start;
	int size = colmode ? nrows : ncols;
	for (i = start; i < size; ++ i) {
		if (colmode) {
			address = i * ncols + pos;
			}
		else {
			address = pos * nrows + i;
			}
		long double value = fabs(pmat[address]);
		if (maxval < value) {
			maxval = value;
			maxpos = i;
			}
		}
	return maxpos;
	}
	
void ldchangerows(long double * pmat, int rows, int ncols, int row1, int row2) {
	int row1pos = row1 * ncols;
	int row2pos = row2 * ncols;
	for (int j = 0; j < ncols; ++ j) {
		long double value = pmat[row1pos + j];
		pmat[row1pos + j] = pmat[row2pos + j];
		pmat[row2pos + j] = value;
		}
	}


// Funções para resover sistemas triangulares em diversas precisões
float * fmtrisolve(float * pmat, int nrows, int ncols, bool superior) {
	float * result;
	result = (float *) malloc(nrows * sizeof(float));
	if (result == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", nrows, ncols);
		return NULL;
		}
	if (superior) {
		for (int i = nrows - 1; i >= 0; -- i) {
			float divisor = pmat[i * (ncols + 1)];
			if (divisor == 0) {
				printf("O sistema é singular! \n");
				return NULL;
				}
			float sum = 0;
			for (int j = i + 1; j < nrows; ++ j) {
				float coef = pmat[i * ncols + j];
				float x = result[j];
				sum += coef * x;
				if (debuglevel_ >= 3) {
					printf("Coluna %d: %f = %f * %f \n", j, sum, coef, x);
					}
				}
			float parm = pmat[(i + 1) * ncols - 1];
			float valor = (parm - sum) / divisor;
			result[i] = valor;
			if (debuglevel_ >= 2) {
				printf("Linha %d: %f = (%f - %f) / %f \n", i, valor, parm, sum, divisor);
				}
			}
		}
	else {
		for (int i = 0; i < nrows; ++ i) {
			float divisor = pmat[i * (ncols + 1)];
			if (divisor == 0) {
				printf("O sistema é singular! \n");
				return NULL;
				}
			float sum = 0;
			for (int j = 0; j < i; ++ j) {
				float coef = pmat[i * ncols + j];
				float x = result[j];
				sum += coef * x;
				if (debuglevel_ >= 3) {
					printf("Coluna %d: %f = %f * %f \n", j, sum, coef, x);
					}
				}
			float parm = pmat[(i + 1) * ncols - 1];
			float valor = (parm - sum) / divisor;
			result[i] = valor;
			if (debuglevel_ >= 2) {
				printf("Linha %d: %f = (%f - %f) / %f \n", i, valor, parm, sum, divisor);
				}
			}
		}
	return result;
	}

double * dmtrisolve(double * pmat, int nrows, int ncols, bool superior) {
	double * result;
	result = (double *) malloc(nrows * sizeof(double));
	if (result == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", nrows, ncols);
		return NULL;
		}
	if (superior) {
		for (int i = nrows - 1; i >= 0; -- i) {
			double divisor = pmat[i * (ncols + 1)];
			if (divisor == 0) {
				printf("O sistema é singular! \n");
				return NULL;
				}
			double sum = 0;
			for (int j = i + 1; j < nrows; ++ j) {
				double coef = pmat[i * ncols + j];
				double x = result[j];
				sum += coef * x;
				if (debuglevel_ >= 3) {
					printf("Coluna %d: %f = %f * %f \n", j, sum, coef, x);
					}
				}
			double parm = pmat[(i + 1) * ncols - 1];
			double valor = (parm - sum) / divisor;
			result[i] = valor;
			if (debuglevel_ >= 2) {
				printf("Linha %d: %f = (%f - %f) / %f \n", i, valor, parm, sum, divisor);
				}
			}
		}
	else {
		for (int i = 0; i < nrows; ++ i) {
			double divisor = pmat[i * (ncols + 1)];
			if (divisor == 0) {
				printf("O sistema é singular! \n");
				return NULL;
				}
			double sum = 0;
			for (int j = 0; j < i; ++ j) {
				double coef = pmat[i * ncols + j];
				double x = result[j];
				sum += coef * x;
				if (debuglevel_ >= 3) {
					printf("Coluna %d: %f = %f * %f \n", j, sum, coef, x);
					}
				}
			double parm = pmat[(i + 1) * ncols - 1];
			double valor = (parm - sum) / divisor;
			result[i] = valor;
			if (debuglevel_ >= 2) {
				printf("Linha %d: %f = (%f - %f) / %f \n", i, valor, parm, sum, divisor);
				}
			}
		}
	return result;
	}
	
long double * ldmtrisolve(long double * pmat, int nrows, int ncols, bool superior) {
	long double * result;
	result = (long double *) malloc(nrows * sizeof(long double));
	if (result == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", nrows, ncols);
		return NULL;
		}
	if (superior) {
		for (int i = nrows - 1; i >= 0; -- i) {
			long double divisor = pmat[i * (ncols + 1)];
			if (divisor == 0) {
				printf("O sistema é singular! \n");
				return NULL;
				}
			long double sum = 0;
			for (int j = i + 1; j < nrows; ++ j) {
				long double coef = pmat[i * ncols + j];
				long double x = result[j];
				sum += coef * x;
				if (debuglevel_ >= 3) {
					printf("Coluna %d: %Lf = %Lf * %Lf \n", j, sum, coef, x);
					}
				}
			long double parm = pmat[(i + 1) * ncols - 1];
			long double valor = (parm - sum) / divisor;
			result[i] = valor;
			if (debuglevel_ >= 2) {
				printf("Linha %d: %Lf = (%Lf - %Lf) / %Lf \n", i, valor, parm, sum, divisor);
				}
			}
		}
	else {
		for (int i = 0; i < nrows; ++ i) {
			long double divisor = pmat[i * (ncols + 1)];
			if (divisor == 0) {
				printf("O sistema é singular! \n");
				return NULL;
				}
			long double sum = 0;
			for (int j = 0; j < i; ++ j) {
				long double coef = pmat[i * ncols + j];
				long double x = result[j];
				sum += coef * x;
				if (debuglevel_ >= 3) {
					printf("Coluna %d: %Lf = %Lf * %Lf \n", j, sum, coef, x);
					}
				}
			long double parm = pmat[(i + 1) * ncols - 1];
			long double valor = (parm - sum) / divisor;
			result[i] = valor;
			if (debuglevel_ >= 2) {
				printf("Linha %d: %Lf = (%Lf - %Lf) / %Lf \n", i, valor, parm, sum, divisor);
				}
			}
		}
	return result;
	}

	
// Funções para cálculo da norma 2 das matrizes em diversas precisões
float fmnorm2(float * pmat, int nrow, int ncol) {
	int size = nrow * ncol;
	float sum = 0;
	while (size -- > 0) {
		float valor = * pmat ++;
		sum += valor * valor;
		}
	return sqrt(sum);
	}

double dmnorm2(double * pmat, int nrow, int ncol) {
	int size = nrow * ncol;
	double sum = 0;
	while (size -- > 0) {
		double valor = * pmat ++;
		sum += valor * valor;
		}
	return sqrt(sum);
	}

long double ldmnorm2(long double * pmat, int nrow, int ncol) {
	int size = nrow * ncol;
	long double sum = 0;
	while (size -- > 0) {
		long double valor = * pmat ++;
		sum += valor * valor;
		}
	return sqrt(sum);
	}

// Funções para leitura das matrizes gravadas pelo MATLAB	
double fgetval(char ** pbuffer) {
// Extrai o primeiro valor existente no "buffer"
	char * pchar = * pbuffer, valor [bufsize], * pvalor = valor;
	while (isblank(* pchar)) {
		pchar ++;
		}
	while (! isblank(* pchar) )	{
		* pvalor ++ = * pchar ++;
		}
	* pvalor = '\0';
	if (debuglevel_ >= 3) {
		printf(" '%s' -> %f", valor, atof(valor));
		}
	* pbuffer = pchar;
	return atof(valor);
	}

double * lermat(const char * fname, int * pnrows, int * pncols) {
// Carrega os dados do arquivo 'fname', gravado pelo MATLAB, numa matriz.
// Retorna a matriz e informa suas dimensões ('nrows' x 'ncols').
// Considera que a matriz está gravada no formato correto.
	// Tenta abrir o arquivo
	char * pchar, * pbuf, buf [bufsize];
	FILE * fp = fopen (fname, "r");
	if (fp == NULL) {
		printf("Não conseguiu ler o arquivo %s! \n", fname);
		return NULL;
		}
	// Despreza as 3 primeiras linhas
	for (int i = 0; i < 3; ++ i) {
		fgets(buf, bufsize, fp);
		}
	// Obtém as dimensões da matriz
	int nrows, ncols;
	fgets(buf, bufsize, fp);
	nrows = atoi(buf + POSROWNBR);
	fgets(buf, bufsize, fp);
	ncols = atoi(buf + POSCOLNBR);
	if (debuglevel_ >= 1) {
		printf("Arquivo %s: linhas = %d, colunas = %d. \n", fname, nrows, ncols);
		}
	double * result, * pval;
	result = pval = (double *) malloc(nrows * ncols * sizeof(double));
	if (pval == NULL) {
		printf("Não conseguiu alocar memória para a matriz %d x %d! \n", nrows, ncols);
		fclose(fp);
		return NULL;
		}
	for (int i = 0; i < nrows; ++ i) {
		buf[bufsize - 2] = '\0';
		fgets(buf, bufsize, fp);
		if (buf[bufsize - 2] != '\0') {
			printf("A matriz tem colunas demais!");
			fclose(fp);
			return NULL;
			}
		pbuf = buf;
		for (int k = 0; k < ncols; ++ k) {
			double valor;
			valor = fgetval(& pbuf);	
			* pval ++ = valor;
			if (debuglevel_ >= 2) {
				printf(" %f ", valor);
				}
			}
		if (debuglevel_ >= 2) {
			printf("\n");
			}
		}
	* pnrows = nrows;
	* pncols = ncols;
	fclose(fp);
	return result;
	}
