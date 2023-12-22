#include <stdlib.h>
#include <math.h>
#include <complex.h>

void Dot(double complex *C, double complex *A, double complex *B, int dim) {
    int i, j, k;
    for(i = 0; i < dim*dim; i++) {
        C[i] = 0;
    }
    for(i = 0; i < dim; i++) {
        for(j = 0; j < dim; j++) {
            for(k = 0; k < dim; k++) {
                C[dim*i + j] += A[dim*i + k] * B[dim*k + j];
            }
        }
    }
}


void Dot_element(double complex *C, double complex *A, double complex *B, int dim) {
    int i;
    for(i = 0; i < dim*dim; i++) {
        C[i] = A[i] * B[i];
    }
}

void Commutator(double complex *C, double complex *A, double complex *B, int dim) {
    int i;
    double complex *AB, *BA;
    AB = (double complex*) calloc(dim*dim, sizeof(double complex));
    BA = (double complex*) calloc(dim*dim, sizeof(double complex));

    Dot(AB, A, B, dim);
    Dot(BA, B, A, dim);
    for(i = 0; i < dim*dim; i++) {
        C[i] = AB[i] - BA[i];
    }
    free(AB);
    free(BA);
}