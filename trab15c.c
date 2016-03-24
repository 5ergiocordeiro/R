#include <stdio.h>
#include <stdlib.h>

int main(int argc, char * argv[]) {
	long min = 1e9, val, dif;
	char mincomb [50] ;
	for (int a=700;a<900;++a) {
		for (int b=700;b<900;++b) {
			for (int c=700;c<900;++c) {
				val = a * c * (2*b - 1);
				dif = val - 1000000000l;
				if ( abs(dif) < min) {
					min = abs(dif);
					sprintf(mincomb,"%d %d %d %d %d", a, b, c, val,dif);
					}
				// printf("%d %d %d -> %d %d\n", a, b, c, val,dif);
				}
			}
		}
	printf("mincomb = %s",mincomb);
	return 0;
	}