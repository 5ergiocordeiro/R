#include <stdio.h>
#include <stdlib.h>

typedef union {
  float f;
  struct {
    unsigned int mantisa : 23;
    unsigned int exponent : 8;
    unsigned int sign : 1;
  } parts;
} double_cast;

int main(int argc, char * argv[]) {
  double_cast d1;
  d1.f = atof(argv[1]);
  printf("sign = %x\n",d1.parts.sign);
  printf("exponent = %x\n",d1.parts.exponent);
  printf("mantisa = %x\n",d1.parts.mantisa);
  return 0;
}