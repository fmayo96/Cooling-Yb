#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "aux.h"

//=================Parameters===============//
//Fundamental constants
const double pi = 3.14159;
const double h = 6.62607015e-34;
const double hbar = h/(2*pi);
const double c = 2.99792458e10;
const double kb = 1.380649e-23;
const double eps0 = 8.8541878128e-14;

//System parameters
const double d = 3.34e-31;
const double spont_em = 314.5;
const double gamma_nr = 1e12;
const double n = 1.45;
const double alpha_imp = 4e-4;
const double alpha_rad = 1e-2;
const int dim = 7;
const double k12 = 6.35e11;
const double k23 = 2.8e11;
const double k34 = 2.25e11;
const double k56 = 5.4e11;
const double k67 = 4.1e11;

// Transition energies
const double w1 = 237*h*c;
const double w2 = 138*h*c;
const double w3 = 102*h*c;
const double w = 9811*h*c;
const double w5 = 132*h*c;
const double w6 = 150*h*c;

//Energy levels
const double E0 = 0;
const double E1 = w1;
const double E2 = E1+w2;
const double E3 = E2+w3;
const double E4 = E3+w;
const double E5 = E4+w5;
const double E6 = E5+w6;
//===============================================//

int main(){
    double complex *state;
    state = (double complex *) calloc(dim*dim, sizeof(double complex));
    


    free(state);
    return 0;
}