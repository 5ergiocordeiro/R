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
#include <time.h>

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
	// Inicializa o gerador de números aleatórios
	time_t t;
	srand((unsigned) time(&t));
	// Inicializa os contadores
	int contador [] = {0, 0, 0, 0, 0, 0, 0, 0};
	// Gera os valores e conta a frequência
	// não é preciso armazenar o valor
	for (int idx = 0; idx < size; ++idx) {
		int val =  rand() % LIMITE_SEQ;
		int pos = val / LIMITE_FAIXA;
		contador[pos]++;
		}
	// Imprime o resultado
	for (int idx = 0; idx < NUM_FAIXAS; ++idx) {
		printf("Faixa %d: %d valores\n",idx+1,contador[idx]);
		}
	return 0;
	}
	
	