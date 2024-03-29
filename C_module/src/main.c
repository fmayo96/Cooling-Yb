#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "aux.h"
#include "time_evol.h"

//=================Parameters===============//
//Fundamental constants
const double pi = 3.14159;
const double h = 6.62607015e-34;
double hbar;
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
double w1;
double w2;
double w3;
double w ;
double w5;
double w6;

//Energy levels
double E0;
double E1;
double E2;
double E3;
double E4;
double E5;
double E6;
//===============================================//

int main(){
    int i;
    double temp = 300, tf =  0.003179650238, dt = 5e-13, j_0_3 = 0, j_0_4 = 0.1, trace, rho_0;
    double complex *state, *H_evo;
    state = (double complex*) calloc(dim*dim, sizeof(double complex));
    H_evo = (double complex*) calloc(dim*dim, sizeof(double complex));

    hbar = h/(2*pi);
    w1 = 237*h*c;
    w2 = 138*h*c;
    w3 = 102*h*c;  
    w = 9811*h*c;
    w5 = 132*h*c;
    w6 = 150*h*c;
    E0 = 0;
    E1 = w1;
    E2 = E1+w2;
    E3 = E2+w3;
    E4 = E3+w;
    E5 = E4+w5;
    E6 = E5+w6;

    double start, end;
	start = clock();

    char filename[255];
    sprintf(filename,"outs/results_temp=%.1lf_j0=%.2lf.txt", temp, j_0_4);

    FILE *fp=fopen(filename,"w");
    
    
    Thermal_state(state, temp);
    Effective_hamiltonian(H_evo, temp, j_0_3, j_0_4);
    Time_evol(state, tf, dt, H_evo);
    
 
    for(i = 0; i < dim-1; i++){
            fprintf(fp, " %e,", creal(state[dim*i+i]));
        }
        fprintf(fp, "%e \n", creal(state[dim*6+6]));
 
    end = clock();
    printf("Program duration: %.2lf \n", (double) (end - start) / CLOCKS_PER_SEC);   
    trace = Trace(state);
    rho_0 = Rho_0(state);
    printf("Trace = %.4lf \n", trace);
    printf("Rho_0 = %.4lf \n", rho_0);
    fclose(fp);
    free(state);
    free(H_evo);
    return 0;
}
