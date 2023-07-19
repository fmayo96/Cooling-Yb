#ifndef POWER_H
#define POWER_H
double Net_power(double temp, double j_0_3, double j_0_4, double tf, double dt);
double Eps_vib_C1(double complex *initial_state, double complex *final_state, double rho_0);
double Eps_vib_C2(double complex *initial_state, double complex *final_state, double rho_0);
double Eps_vib_H1(double complex *initial_state, double complex *final_state, double rho_0);
double Eps_vib_H2(double eps_vib_h1, double eps_vib_c1);
#endif