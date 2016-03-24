/*
path=%path%;C:\Octave\Octave-4.0.0\bin
Solução de sistema de equações lineares (Ax = b) por meio dos solvers interativos da biblioteca Sparskit.
Uso:
	solve Afile bfile maxiter erro prec solver
onde
	Afile é o nome do arquivo contendo os valores da matriz A
	bfile é o nome do arquivo contendo os valores do vetor b
	maxiter é o número máximo de iterações permitido
	erro é o erro absoluto tolerado
	prec é o nome do precondicionador a usar
	solver é o nome do solver a ser usado
Para mais informações, consultar o código fonte da biblioteca (arquivos iters.f, ilut.f, matvec.f e formats.f e blassm.f).
Listas dos solvers e precondicionadores permitidos encontram-se no código abaixo.
Os arquivos formats.f e blassm.f contêm muitas dependências, por isso apenas a parte que interessava para esta aplicação foi mantida.
Limitações:
	para o número de condição foi empregada uma fórmula aproximada,
	não foram utilizadas todas as funcionalidades dos precondicionadores
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define bufsize 5000
#define size 37
#define numsolvers 10
#define numprecs 8
#define Krylov 15
// #define TESTA_PREC


// Assinatura redundante apenas para compatibilidade com a biblioteca Sparskit
extern"C" { double distdot_(int *, double *, int *, double *, int *); }
extern"C" { double ddot_(int *, double *, int *, double *, int *); }
extern"C" { void daxpy_(int *, double *, double *, int *, double *, int *); }
extern"C" { double dnrm2_(int *, double *, int *); }

// Funções da biblioteca Sparskit (em Fortran)
// ... solvers
extern"C" { void cg_(int *, double *, double *, int *, double *, double *); }
extern"C" { void cgnr_(int *, double *, double *, int *, double *, double *); }
extern"C" { void bcg_(int *, double *, double *, int *, double *, double *); }
extern"C" { void dbcg_(int *, double *, double *, int *, double *, double *); }
extern"C" { void bcgstab_(int *, double *, double *, int *, double *, double *); }
extern"C" { void tfqmr_(int *, double *, double *, int *, double *, double *); }
extern"C" { void fom_(int *, double *, double *, int *, double *, double *); }
extern"C" { void gmres_(int *, double *, double *, int *, double *, double *); }
extern"C" { void fgmres_(int *, double *, double *, int *, double *, double *); }
extern"C" { void dqgmres_(int *, double *, double *, int *, double *, double *); }
// ... BLAS
extern"C" { void amux_(int *, double *, double *, double *, int *, int *); }
extern"C" { void atmux_(int *, double *, double *, double *, int *, int *); }
extern"C" { void amub_(int *, int *, int *, double *, int *, int *, double *, int *, int *, double *, int *, int *, int *, int *, int *); }
// ... conversões
extern"C" { void dnscsr_(int *, int *, int *, double *, int *, double *, int *, int *, int *); }
extern"C" { void csrdns_(int *, int *, double *, int *, int *, double *, int *, int *); }
extern"C" { void csrmsr_(int *, double *, int *, int *, double *, int *, double *, int *); }
extern"C" { void msrcsr_(int *, double *, int *, double *, int *, int *, double *, int *); }
// ... ILU
extern"C" { void ilut_(int *, double *, int *, int *, int *, double *, double *, int *, int *, int *, double *, int *, int *); }
extern"C" { void ilud_(int *, double *, int *, int *, double *, double *, double *, int *, int *, int *, double *, int *, int *); }
extern"C" { void ilutp_(int *, double *, int *, int *, int *, double *, double *, int *, double *, int *, int *, int *, double *, int *, int *, int *); }
extern"C" { void iludp_(int *, double *, int *, int *, double *, double *, double *, int *, double *, int *, int *, int *, double *, int *, int *, int *); }
extern"C" { void iluk_(int *, double *, int *, int *, int *, double *, int *, int *, int *, int *, double *, int *, int *); }
extern"C" { void ilu0_(int *, double *, int *, int *, double *, int *, int *, int *, int *); }
extern"C" { void milu0_(int *, double *, int *, int *, double *, int *, int *, int *, int *); }
extern"C" { void lusol_(int *, double *, double *, double *, int *, int *); }
extern"C" { void lutsol_(int *, double *, double *, double *, int *, int *); }
extern"C" { void pgmres_(int *, int *, double *, double *, double *, double *, int *, int *, double *, int *, int *, double *, int *, int *, int *); }

// Wrappers para chamada de precondicionadores da biblioteca
void myilut(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void myilud(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void myilutp(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void myiludp(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void myiluk(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void myilu0(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
void mymilu0(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);

// Cálculo do espaço de Krylov
int Ksize1();
int Ksize2();
int Ksize3();


typedef void SolverFn(int *, double *, double *, int *, double *, double *);
typedef int KrylovFn(void);
typedef void PrecFn(int *, double *, int *, int *, double *, int *, int *, int *, double *, int *, int *);
typedef struct {
	const char * nome;
	SolverFn * function;
	int wsize;
	KrylovFn * kfn ;
	} SolverInfo;
typedef struct {
	const char * nome;
	PrecFn * function;
	int mode;
	} PrecInfo;


// Solvers implementados	
static SolverInfo Solver[] = {
	{ "CG", 		cg_,		5, 	NULL, },
	{ "CGNR",		cgnr_,		5, 	NULL, },
	{ "BCG",		bcg_,		7, 	NULL, },
	{ "DBCG",		dbcg_,		11, NULL, },
	{ "BCGSTAB",	bcgstab_,	8, 	NULL, },
	{ "TFQMR",		tfqmr_,		11, NULL, },
	{ "FOM",		fom_,		0, 	Ksize1, },
	{ "GMRES",		gmres_,		0, 	Ksize1, },
	{ "FGMRES",		fgmres_,	0, 	Ksize2, },
	{ "DQGMRES",	dqgmres_,	0, 	Ksize3, }
	};

// Precondicionadores implementados
static PrecInfo Preconditioner[] = {
	{ "0", 			NULL,		0, },
	{ "ILUT",		myilut,		2, },
	{ "ILUTP",		myilutp,	2, },
	{ "ILUD",		myilud,		2, },
	{ "ILUDP",		myiludp,	2, },
	{ "ILUK",		myiluk,		2, },
	{ "ILU0",		myilu0,		2, },
	{ "MILU0",		mymilu0,	2, }
	};

	
// Funções de apoio (BLAS nível 1)
double distdot_(int * psize, double * x, int * ix, double * y, int * iy) {
// Calcula o produto escalar de dois vetores.
// Assinatura pre-definida pela biblioteca Sparskit.
	double sum = 0;
	int vsize = * psize;
	if (* ix != 1 || * iy != 1) {
		printf("Situação não tratada (inesperada) no cálculo do produto escalar");
		exit(1);
		}
	while (vsize -- > 0) {
		sum += (* x ++) * (* y ++);
		}
	return sum;
	}

double ddot_(int * psize, double * x, int * ix, double * y, int * iy) {
// Calcula o produto escalar de dois vetores.
// Assinatura pre-definida pela biblioteca Sparskit.
	return distdot_(psize, x, ix, y, iy) ;
	}

double dnrm2_(int * psize, double * x, int * ix) {
// Calcula a norma Euclideana de um vetor
// Assinatura pre-definida pela biblioteca Sparskit.
    double sum = 0;
	int vsize = * psize;
	if (* ix != 1) {
		printf("Situação não tratada (inesperada) no cálculo do produto escalar");
		exit(1);
		}
	while (vsize -- > 0) {
		sum += (* x) * (* x);
		++ x;
		}
	return sqrt(sum);
	}

void daxpy_(int * psize, double * pa, double * x, int * ix, double * y, int * iy) {
// Calcula y = ax + y, onde a é um escalar e x e y são vetores
// Assinatura pre-definida pela biblioteca Sparskit.
    double sum = 0;
	int vsize = * psize, a = * pa;
	if (* ix != 1 || * iy != 1) {
		printf("Situação não tratada (inesperada) no cálculo do produto escalar");
		exit(1);
		}
	while (vsize -- > 0) {
		* y ++ += a * (* x ++);
		}
	}

	
// Funções para leitura de dados em disco
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

	
// Funções específicas
// Wrappers para chamada de precondicionadores da biblioteca
void myilut(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
	int lfil = 5;
	double droptol = 1e-3;
    ilut_(n, a, ja, ia, & lfil, & droptol, alu, jlu, ju, iwk, w, jw, ierr);
	}

void myilud(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
	double alph = 0.5, tol = 1e-3;
    ilud_(n, a, ja, ia, & alph, & tol, alu, jlu, ju, iwk, w, jw, ierr);
	}

void myilutp(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
	int mbloc = * n, lfil = 5;
	int * iperm = (int *) malloc(2 * (*n) * sizeof(int));
	if (iperm == NULL) {
		printf("Não conseguiu alocar a memória de trabalho necessária para o precondicionador");
		exit(1);
		}
	double permtol = 0, droptol = 1e-3;
    ilutp_(n, a, ja, ia, & lfil, & droptol, & permtol, & mbloc, alu, jlu, ju, iwk, w, jw, iperm, ierr);
	}

void myiludp(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
	int mbloc = * n;
	int * iperm = (int *) malloc(2 * (*n) * sizeof(int));
	if (iperm == NULL) {
		printf("Não conseguiu alocar a memória de trabalho necessária para o precondicionador");
		exit(1);
		}
	double alph = 0.5, tol = 1e-3, permtol = 0;
    iludp_(n, a, ja, ia, & alph, & tol, & permtol, & mbloc, alu, jlu, ju, iwk, w, jw, iperm, ierr);
	}

void myiluk(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
	int lfil = 5, * levs = (int *) malloc((*iwk) * sizeof(int));
	if (levs == NULL) {
		printf("Não conseguiu alocar a memória de trabalho necessária para o precondicionador");
		exit(1);
		}
    iluk_(n, a, ja, ia, & lfil, alu, jlu, ju, levs, iwk, w, jw, ierr);
	}

void myilu0(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
    ilu0_(n, a, ja, ia, alu, jlu, ju, jw, ierr);
	}

void mymilu0(int * n, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju, int * iwk, double * w, int * jw, int * ierr) {
    milu0_(n, a, ja, ia, alu, jlu, ju, jw, ierr);
	}


void teste_ilu(double * b, double * x, double maxerr, int maxiter, double * a, int * ja, int * ia, double * alu, int * jlu, int * ju) {
// Testa o resultado do precondicionamento.
// Apenas para fins de desenvolvimento.
	int zero = 1, n = size, retcode, wksize = Krylov;
	double * pvv = (double *) malloc((wksize + 1) * size * sizeof(double));
	if (pvv == NULL) {
		printf("Não conseguiu alocar a memória de trabalho necessária para o solver de teste");
		exit(1);
		}
	pgmres_(& n, & wksize, b, x, pvv, & maxerr, & maxiter, & zero, a, ja, ia, alu, jlu, ju, & retcode);
	if (retcode != 0) {
		printf("Erro no cálculo do solver de teste PGMRES: %d", retcode);
		exit(1);
		}
	printf("Solução do solver de teste PGMRES: [");
	for (int i = 0; i < size ; ++ i) {
		if (i > 0) {
			printf(", ");
			}
		printf("%f", x[i]);
		}
	}


void getprec(double * a, int * ja, int * ia, PrecInfo * pprec, double * alu, int * jlu, int * ju) {
// Obtém o precondicionador em formato MSR (modified sparse row)
	double * wp = (double *) malloc((size + 1) * sizeof(double));
	int * jwp = (int *) malloc(2 * size * sizeof(int));
	int n = size, n2 = n * n + 100, retcode;
	(pprec -> function) (& n, a, ja, ia, alu, jlu, ju, & n2, wp, jwp, & retcode);
	if (retcode != 0) {
		printf("Erro no cálculo do precondicionador: %d", retcode);
		exit(1);
		}
	free(wp); free(jwp);
	}

	
// Cálculo do espaço de Krylov
int Ksize1() {
	return (size + 3) * (Krylov + 2) + Krylov * ( Krylov + 1) / 2.0 + 1;
	}

int Ksize2() {
	return 2 * size * (Krylov + 1) + Krylov * ( Krylov + 1) / 2.0 + 3 * Krylov + 3;
	}

int Ksize3() {
	return size + (Krylov + 1) * (2 * size + 4);
	}

int solve(double * a, int * ja, int * ia, double * b, double * x, int maxiter, double maxerr, int idxprec, int idxsolver) {
// Resolve o sistema linear ax = b
// a está na forma CSR
// Utiliza diversas combinações de solvers e precondicionadores
	// Inicializa o precondicionador, se houver
	PrecInfo * pprec = & Preconditioner[idxprec];	
	double * alu;
	int * jlu, * ju;
	if (idxprec > 0) {
		alu = (double *) malloc(size * size * sizeof(double));
		jlu = (int *) malloc(size * size * sizeof(int));
		ju = (int *) malloc(size * sizeof(int));
		if (alu == NULL || jlu == NULL || ju == NULL) {
			printf("Não conseguiu alocar a memória de trabalho necessária para o precondicionador");
			exit(1);
			}
		getprec(a, ja, ia, pprec, alu, jlu, ju);
#ifdef TESTA_PREC
		teste_ilu(b, x, maxerr, maxiter, a, ja, ia, alu, jlu, ju);
		exit(0);
#endif
		}
	// Inicializa o solver
	int retcode, n = size, one = 1;
	int ipar[16];
	double fpar[16], flop;
	SolverInfo * psolver = & Solver[idxsolver];
	int wksize = size * psolver -> wsize;
	if (wksize <= 0) {
		wksize = (psolver -> kfn)();
		}
	double * ws = (double *) malloc(wksize * sizeof(double));
	if (ws ==  NULL) {
		printf("Não conseguiu alocar a memória de trabalho necessária para o solver");
		return 1;
		}
	fpar[0] = maxerr;
	fpar[1] = maxerr * 1e-4;
	fpar[10] = flop = 0;
	ipar[0] = 0;
	ipar[1] = pprec -> mode;
	ipar[2] = 1;
	ipar[3] = wksize;
	ipar[4] = Krylov;
	ipar[5] = maxiter;
	ipar[6] = 0;
	// Executa o solver interativamente
	do {
		(psolver -> function) (& n, b, x, ipar, fpar, ws);
		retcode = ipar[0];
		switch (retcode) {
			case 1:
				amux_(& n, ws + ipar[7] - 1, ws + ipar[8] - 1, a, ja, ia ) ;
				flop += (2 * size - 1) * size;
				printf(".");
				break;
			case 2:
				atmux_(& n, ws + ipar[7] - 1, ws + ipar[8] - 1, a, ja, ia ) ;
				flop += (2 * size - 1) * size;
				printf("|");
				break;			
			case 3:
			case 5:
				if (idxprec == 0) {
					printf("Situação não tratada (inesperada): condicionador chamado");
					exit(1);
					}
				lusol_(& n, ws + ipar[7] - 1, ws + ipar[8] - 1, alu, jlu, ju);
				flop += (2 * size - 1) * size;
				printf("+");
				break;
			case 4:
			case 6:
				if (idxprec == 0) {
					printf("Situação não tratada (inesperada): condicionador chamado");
					exit(1);
					}
				lutsol_(& n, ws + ipar[7] - 1, ws + ipar[8] - 1, alu, jlu, ju);
				flop += (2 * size - 1) * size;
				printf("-");
				break;	
			default:
				if (retcode > 0) {
					printf("Situação não tratada (inesperada): retcode = %d\n", retcode);
					exit(1);
					}
				}
		} while (retcode > 0);
	if (retcode < 0) {
		printf("Não obteve a solução após %d iterações. Erro = %d\n", ipar[6], retcode);
		return 1;
		}
	// Exibe a solução
	double ka, val, max = 0, min = 1e9;
	double gflop = (flop + fpar[10]) / 1e9;
	printf("Solução após %d iterações (%f GFlops): [", ipar[6], gflop);
	for (int i = 0; i < size ; ++ i) {
		if (i > 0) {
			printf(", ");
			}
		printf("%f", x[i]);
		}
	if (idxprec == 0) {
		ka = max / min;
		}
	if (idxprec == 1) {
		ka = 1;
		}
	printf("] \nErro < %f, ka = %f", fpar[5], ka);
	return 0;
	}


int findsolver(const char * nome) {
	for (int i = 0; i < numsolvers; ++ i) {
		if (! strcmp (nome, Solver[i].nome))
			return i;
		}
	return -1;
	}


int findprec(const char * nome) {
	for (int i = 0; i < numprecs; ++ i) {
		if (! strcmp (nome, Preconditioner[i].nome))
			return i;
		}
	return -1;
	}
	
	
int convert(double * A, double ** pA, int ** pja, int ** pia) {
// Converte uma matriz para o formato CSR (compressed sparse row)
	int retcode, n = size, nmax = 0;
	* pA = A + size * size - 1;
	while (* pA >= A) {
		nmax += (** pA != 0) ;
		-- (* pA);
		}
	* pA = (double *) malloc(nmax * sizeof(double));
	* pja = (int *) malloc(nmax * sizeof(int));
	* pia = (int *) malloc((size + 1) * sizeof(int));	
	if (* pA == NULL || * pja == NULL || * pia == NULL) {
		printf("Não conseguiu alocar memória para a matriz convertida");
		exit(1);
		}		
	dnscsr_(& n, & n, & nmax, A, & n, * pA, * pja, * pia, & retcode);
	if (retcode != 0) {
		printf("Erro %d ao converter matriz para o formato CSR", retcode);
		exit(1);
		}
	return nmax;
	}

	
int main(int argc, char * argv[]) {
// Solução de sistema de equações lineares (Ax = b) por meio dos solvers da biblioteca Sparskit. 
	// Carrega os dados a partir dos arquivos em disco
	double A[size][size], b[size], x[size];		// Ax = b
	int noread;
	noread = load(argv[1], &A[0][0], size, size);
	if (noread > 0) {
		return 1;
		}
	noread = load(argv[2], b, size, 1);
	if (noread > 0) {
		return 1;
		}
	// Lê os demais parâmetros
	int maxiter = atoi(argv[3]) ;
	if (maxiter <= 0) {
		printf("Número máximo de iterações inválido: %d", maxiter);
		return 1;
		}
	double erro = atof(argv[4]);
	if (erro <= 0) {
		printf("Erro tolerado inválido: %f", erro);
		return 1;
		}
	int idxprec = findprec(argv[5]);
	if (idxprec < 0) {
		printf("Precondicionador desconhecido: %s", argv[5]);
		return 1;
		}
	int idxsolver = findsolver(argv[6]);
	if (idxsolver < 0 ) {
		printf("Solver desconhecido: %s", argv[6]);
		return 1;
		}
	// Inicializa o vetor solução (x0)
	for (int i = 0; i < size ; ++ i) {
		// (teoricamente, pode ser qualquer coisa)
		x[i] = 0;
		}
	// Converte a matriz de coeficientes para o formato CSR
	double * pA;
	int * pja, * pia;
	int nmax = convert(& A[0][0], & pA, & pja, & pia);
	printf("nmax = %d (%s esparsa)\n", nmax, (nmax * 1.0 / (size * size)) > 0.5 ? "pouco" : "muito");
	// Resolve o sistema
	solve(pA, pja, pia, b, x, maxiter, erro, idxprec, idxsolver);
	}