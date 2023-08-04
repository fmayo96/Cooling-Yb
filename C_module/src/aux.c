#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"


void Thermal_state(double complex *thermal_state, double temp){
    int i;
    double beta = (double)1/(kb*temp), Z = 0;
    *(thermal_state + dim*0 + 0) = exp(-beta*E0);
    *(thermal_state + dim*1 + 1) = exp(-beta*E1);
    *(thermal_state + dim*2 + 2) = exp(-beta*E2);
    *(thermal_state + dim*3 + 3) = exp(-beta*E3);
    *(thermal_state + dim*4 + 4) = exp(-beta*E4);
    *(thermal_state + dim*5 + 5) = exp(-beta*E5);
    *(thermal_state + dim*6 + 6) = exp(-beta*E6);
    for(i=0; i<dim; i++){
        Z += *(thermal_state + dim*i + i);
    }
    for(i=0; i<dim; i++){
        *(thermal_state + dim*i + i) /= Z;
    }
}

double Thermal_num(double temp, double E){
    double beta = (double)1/(kb*temp);
    double thermal_num;
    thermal_num = (double)1/(exp(beta*E) - 1);
    return thermal_num; 
}

double Rabifreq(double j_0){
    double E_0;
    j_0 *= 1e6;
    E_0 = sqrt((double)2*j_0/(c*n*eps0));
    return (double) d*E_0/hbar;
}

double Trace(double complex *state){
    int i;
    double trace = 0;
    for(i = 0; i < dim; i++){
        trace += creal(*(state + dim*i + i));
    }
    return trace;
}

double Rho_0(double complex *state){
    int i;
    double rho_0 = 0;
    for(i = 4; i < dim; i++){
        rho_0 += creal(*(state + dim*i + i));
    }
    return rho_0;
}
