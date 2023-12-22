#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATION_H

void Dot(double complex *C, double complex *A, double complex *B, int dim);
void Dot_element(double complex *C, double complex *A, double complex *B, int dim);
void Commutator(double complex *C, double complex *A, double complex *B, int dim);


#endif