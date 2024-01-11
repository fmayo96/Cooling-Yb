#ifndef TIME_EVOL_H
#define TIME_EVOL_H
void Time_evol(double complex *state, double tf, double dt, double temp, double j_0_3, double j_0_4);
/* void RK4_step(double complex *state, double dt, double temp, double j_0_3, double j_0_4);*/
void RK2_step(double complex *state, double dt, double temp, double j_0_3, double j_0_4);
void Diff(double complex *diff_state, double complex *state, double temp, double j_0_3, double j_0_4);
void Unit_evo(double complex *unit_state, double complex *state, double temp, double j_0_3, double j_0_4);
void Decoh(double complex *decoh_state, double complex *state);
void NonRad(double complex *non_rad_state, double complex *state);
void SpontEm(double complex *spont_em_state, double complex *state);


#endif
