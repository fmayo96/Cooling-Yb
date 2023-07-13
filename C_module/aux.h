#ifndef AUX_H
#define AUX_H
void Thermal_state(double complex *thermal_state, double temp);
double Thermal_num(double temp, double E);
double Rabifreq(double j_0);
double Trace(double complex *state);
#endif