#ifndef TIME_EVOL_H
#define TIME_EVOL_H
void Time_evol(double complex *state, double tf, double dt, double temp, double j_0_3, double j_0_4);
void RK4_step(double complex *state, double dt, double temp, double j_0_3, double j_0_4);
void Diff(double complex *diff_state, double complex*state, double temp, double j_0_3, double j_0_4);
#endif