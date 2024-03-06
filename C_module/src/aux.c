#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"


void Thermal_state(double complex *thermal_state, double temp) {
    int i;
    double beta = (double)1/(kb*temp), Z = 0;
    *(thermal_state + dim*0 + 0) = exp(-beta*E0);
    *(thermal_state + dim*1 + 1) = exp(-beta*E1);
    *(thermal_state + dim*2 + 2) = exp(-beta*E2);
    *(thermal_state + dim*3 + 3) = exp(-beta*E3);
    *(thermal_state + dim*4 + 4) = 0;
    *(thermal_state + dim*5 + 5) = 0;
    *(thermal_state + dim*6 + 6) = 0;
    for(i=0; i<dim; i++){
        Z += *(thermal_state + dim*i + i);
    }
    for(i=0; i<dim; i++){
        *(thermal_state + dim*i + i) /= Z;
    }
}

double Thermal_num(double temp, double E) {
    double beta = (double)1/(kb*temp);
    double thermal_num;
    thermal_num = (double)1/(exp(beta*E) - 1);
    return thermal_num; 
}

double Rabifreq(double j_0) {
    double E_0;
    j_0 *= 1e6;
    E_0 = sqrt((double)2*j_0/(c*n*eps0));
    return (double) d*E_0/hbar;
}

double Trace(double complex *state) {
    int i;
    double trace = 0;
    for(i = 0; i < dim; i++){
        trace += creal(*(state + dim*i + i));
    }
    return trace;
}

double Rho_0(double complex *state) {
    int i;
    double rho_0 = 0;
    for(i = 4; i < dim; i++){
        rho_0 += creal(*(state + dim*i + i));
    }
    return rho_0;
}

void Effective_hamiltonian(double complex *H_evo, double temp, double j_0_3, double j_0_4) {
    int i;
    double Rabi4, Rabi3;
    double *ks, *ws;
    ks = (double*) calloc(6, sizeof(double));
    ws = (double*) calloc(6, sizeof(double));
    Rabi4 = Rabifreq(j_0_4);
    Rabi3 = Rabifreq(j_0_3);
    ks[0] = k12; ks[1] = k23; ks[2] = k34; ks[3] = Rabi4; ks[4] = k56; ks[5] = k67;  
    ws[0] = w1; ws[1] = w2; ws[2] = w3; ws[3] = w; ws[4] = w5; ws[5] = w6; 
    
    for(i = 0; i < dim - 1; i++) {
        H_evo[dim*i + i+1] = ks[i] * sqrt(Thermal_num(temp, ws[i]));
    }
    for(i = 1; i < dim; i++) {
        H_evo[dim*i + i - 1] = ks[i-1] * sqrt(Thermal_num(temp, ws[i-1]));
    }
    H_evo[dim*3 + 4] = Rabi4;
    H_evo[dim*4 + 3] = Rabi4;
    H_evo[dim*4 + 2] = Rabi3;
    H_evo[dim*2 + 2] = Rabi3;
    free(ks);
    free(ws);
}