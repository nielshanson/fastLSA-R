#include <omp.h>
#include <R.h>

// arguments
// m: a square matrix
// n: number of rows/columns of m
// k: the subdiagonal index 
// result: space for requested sub-diagonal

extern "C" void subdiag(double *m, int *n, int *k, double *result) {
	
	int nval = *n, kval = *k;
	int stride = nval + 1;
	for (int i = 0, j = kval; i < nval-kval; ++i, j+= stride) {
		result[i] = m[j];
		
	}
}