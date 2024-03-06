#ifndef AUX_H
#define AUX_H
void Thermal_state(double complex *thermal_state, double temp);
double Thermal_num(double temp, double E);
double Rabifreq(double j_0);
double Trace(double complex *state);
double Rho_0(double complex *state);
void Effective_hamiltonian(double complex *H_evo, double temp, double j_0_3, double j_0_4);
#endif