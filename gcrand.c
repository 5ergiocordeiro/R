/*
gcrand.c
Gera sequência de números aleatórios e calcula a frequência absoluta de incidência em cada faixa de valores.
Uso:
	gcrand size
	onde size é o tamanho da sequência
Testado em GNU C sobre Linux (Ubuntu 12.04.5).
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <windows.h>

#define LIMITE_SEQ		255
#define LIMITE_FAIXA	32
#define NUM_FAIXAS		8

int main(int argc, char * argv[]) {
	// Obtém o tamanho da sequência
	int size = atoi(argv[1]);
	if ( size <= 0 ) {
		printf("Tamanho da sequencia deve ser positivo!\n");
		return 1;
		}
	struct timeval t1, t2;
    HANDLE hProcess = GetCurrentProcess();
    FILETIME ftCreation, ftExit, ftKernel, ftUser1, ftUser2;
    SYSTEMTIME stUser1,stUser2;
    gettimeofday(&t1, NULL);
    GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser1);
	// Inicializa o gerador de números aleatórios
	time_t t;
	srand((unsigned) time(&t));
	// Inicializa os contadores
	int contador [] = {0, 0, 0, 0, 0, 0, 0, 0};
	// Gera os valores e conta a frequência
	// não é preciso armazenar o valor
	for (int i=0;i<1000;++i) {
	for (int idx = 0; idx < size; ++idx) {
		int val =  rand() % LIMITE_SEQ;
		int pos = val / LIMITE_FAIXA;
		contador[pos]++;
		}
	}
	// Imprime o resultado
	
	for (int idx = 0; idx < NUM_FAIXAS; ++idx) {
		printf("Faixa %d: %d valores\n",idx+1,contador[idx]);
		}
	//
    gettimeofday(& t2, NULL);
    GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser2);
    FileTimeToSystemTime(& ftUser1, & stUser1);
    FileTimeToSystemTime(& ftUser2, & stUser2);
    printf("Time = %d.%d, User time = %d usecs \n", 
		t2.tv_sec - t1.tv_sec,t2.tv_usec - t1.tv_usec,
		stUser2.wMilliseconds - stUser1.wMilliseconds);
	return 0;
	}
	
	