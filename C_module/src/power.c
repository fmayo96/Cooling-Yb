#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "constants.h"
#include "aux.h"
#include "time_evol.h"
#include "power.h"


double Net_power(double temp, double j_0_3, double j_0_4, double tf, double dt){
    double complex *initial_state, *final_state;
    double rho_0, n_ion = 1.3e20, N_e = 13, eta_e = 0.92, W_nr = 1.45;
    double eps_vib_c1, eps_vib_c2, eps_vib_h1, eps_vib_h2;
    double pow_cool, pow_nr, pow_imp, pow_net;
    initial_state = (double complex *) calloc(dim*dim, sizeof(double complex));
    final_state = (double complex *) calloc(dim*dim, sizeof(double complex));
    Thermal_state(initial_state, temp);    
    Thermal_state(final_state, temp);  
    Time_evol(final_state, tf, dt, temp, j_0_3, j_0_4);  
    rho_0 = Rho_0(final_state);
    printf("%e \n",rho_0);
    eps_vib_c1 = Eps_vib_C1(initial_state, final_state, rho_0);
    eps_vib_c2 = Eps_vib_C2(initial_state, final_state, rho_0);
    eps_vib_h1 = Eps_vib_H1(initial_state, final_state, rho_0);
    eps_vib_h2 = Eps_vib_H2(eps_vib_h1, eps_vib_c1);

    pow_cool = eta_e * n_ion * N_e * rho_0 * spont_em * (eps_vib_c1 + eps_vib_c2);
    pow_nr = n_ion * N_e * rho_0 * (spont_em * eps_vib_h1 + W_nr * eps_vib_h2);
    pow_imp = (j_0_3 + j_0_4) * alpha_imp * 1e6;
    pow_net = pow_nr + pow_imp - pow_cool;

    printf("Pow_nr = %lf \n", pow_nr);
    printf("Pow_imp = %lf \n", pow_imp);
    printf("Pow_cool = %lf \n", pow_cool);


    free(initial_state);
    free(final_state);
    
    return pow_net;
}
double Eps_vib_C1(double complex *initial_state, double complex *final_state, double rho_0){
    double eps_vib_c1 = E3 + (double)(E0 * (*final_state - *initial_state) + E1 * (*(final_state + dim*1 + 1) - *(initial_state + dim*1 + 1)) +
    E2 * (*(final_state + dim*2 + 2) - *(initial_state + dim*2 + 2)) + E3 * (*(final_state + dim*3 + 3) - *(initial_state + dim*3 + 3))) / rho_0;

    return eps_vib_c1;
}
double Eps_vib_C2(double complex *initial_state, double complex *final_state, double rho_0){
    double eps_vib_c2 = (double)(E4 * (*(final_state + dim*4 + 4) - *(initial_state + dim*4 + 4)) +
    E5 * (*(final_state + dim*5 + 5) - *(initial_state + dim*5 + 5)) + E6 * (*(final_state + dim*6 + 6) - *(initial_state + dim*6 + 6))) / rho_0
    + E3 - E4 - 0.25*(E0 + E1 + E2 + E3);
    return eps_vib_c2;
}
double Eps_vib_H1(double complex *initial_state, double complex *final_state, double rho_0){
    double eps_vib_h1 = E3 + (double)(E0 * (*final_state - *initial_state + 0.25*rho_0) + E1 * (*(final_state + dim*1 + 1) - *(initial_state + dim*1 + 1) + 0.25*rho_0) +
    E2 * (*(final_state + dim*2 + 2) - *(initial_state + dim*2 + 2) + 0.25*rho_0) + E3 * (*(final_state + dim*3 + 3) - *(initial_state + dim*3 + 3) + 0.25*rho_0)) / rho_0;
    return eps_vib_h1;
}
double Eps_vib_H2(double eps_vib_h1, double eps_vib_c1){
    double eps_vib_h2 = E4 - E3 + eps_vib_h1 + eps_vib_c1;
    return eps_vib_h2;
}