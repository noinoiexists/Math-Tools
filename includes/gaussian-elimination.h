#include <inttypes.h>

typedef uint8_t Bit;

double *GaussianElimination(size_t n, double A[n][n+1]);
Bit *GF2GaussianElimination(size_t n, Bit A[n][n+1]);
