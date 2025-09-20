#include <stdio.h>
#include <limits.h>
#include "gaussian-elimination.h"

void swap_rows(size_t n, double A[n][n+1], size_t row1, size_t row2);

double *GaussianElimination(size_t n, double A[n][n+1]){
	size_t min_row;
	for (size_t i=0; i<n; ++i){
		// Pivoting
		min_row=i;
		for (size_t j=0; i+j<n; ++j){
			if (MIN(ABS(A[min_row][i]),ABS(A[i+j][i])) != 0.0 && MIN(ABS(A[min_row][i]),ABS(A[i+j][i])) != ABS(A[min_row][i]))
				min_row = i+j;
		}
		if (ABS(A[min_row][i]) == 0.0){
			fprintf(stderr, "No or Infinite solution\n");
			exit(1);
		} 
		
		if (i != min_row)
			swap_rows(n, A, i, min_row);
		// Forward Elimination
		for (size_t j=1; i+j<n; ++j){

			double factor = A[i+j][i] / A[i][i];
			for (size_t k=0; i+k<=n; ++k){
				A[i+j][i+k] -= factor*A[i][i+k];
			}
		}

	}
	// Back Substitution
	double *soln = malloc(sizeof(double)*n);
	if (!soln){
		fprintf(stderr, "Memory allocation failed\n");
		exit(1);
	}
	for (int32_t i=n-1; i>=0; --i){
		double C = 0.0;
		for (size_t j=1; i+j<n; ++j){
			C += A[i][i+j]*soln[i+j];
		}
		soln[i] = (A[i][n] - C) / A[i][i];
	}
	
	return soln;
}

void swap_rows(size_t n, double A[n][n+1], size_t row1, size_t row2) {
	for(size_t i=0; i<=n; ++i){
		double temp = A[row1][i];
		A[row1][i] = A[row2][i];
		A[row2][i] = temp;
	}
}