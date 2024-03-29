#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "aux.h"
#include "matrix_operations.h"
#include "time_evol.h"

void Time_evol(double complex *state, double tf, double dt, double complex *H_evo){
    long N = (long) (tf/dt), i;
    for(i = 0; i < N; i++){
        Euler_step(state, dt, H_evo);
    }
}

void Euler_step(double complex *state, double dt, double complex *H_evo) {
    int i,j;
    double complex *diff_state;
    diff_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    
    Diff(diff_state, state, H_evo);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            state[dim*i + j] += diff_state[dim*i + j] * dt;
        }
    }
    free(diff_state);
}

void Diff(double complex *diff_state, double complex *state, double complex *H_evo){
    int i;
    double complex *unit_state, *decoh_state, *non_rad_state, *spont_em_state;
    unit_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    decoh_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    non_rad_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    spont_em_state = (double complex*) calloc(dim*dim, sizeof(double complex));

    Unit_evo(unit_state, state, H_evo);
    Decoh(decoh_state, state);
    NonRad(non_rad_state, state);
    SpontEm(spont_em_state, state);

    for(i = 0; i < dim*dim; i++) {
        diff_state[i] = unit_state[i] + decoh_state[i] + non_rad_state[i] + spont_em_state[i];
    }

    free(unit_state);
    free(decoh_state);
    free(non_rad_state);
    free(spont_em_state);
}


void Unit_evo(double complex *unit_state, double complex *state, double complex *H_evo) {
    int i;

    Commutator(unit_state, H_evo, state, dim);
    for(i = 0; i < dim*dim; i++) {
        unit_state[i] *= -I;
    }
}

void Decoh(double complex *decoh_state, double complex *state) {
    int i, j;
    for(i = 0; i < 4; i++) {
        for(j = 4; j < dim; j++) {
            decoh_state[dim*i + j] = -0.5 * gamma_nr * state[dim*i + j]; 
            decoh_state[dim*j + i] = -0.5 * gamma_nr * state[dim*j + i]; 
        }
    }
    for(i = 0; i < dim-1; i++) {
        decoh_state[dim*i + i+1] = -0.5 * gamma_nr * state[dim*i + i+1];
        decoh_state[dim*(i+1) + i] = -0.5 * gamma_nr * state[dim*(i+1) + i];
    } 
}

void NonRad(double complex *non_rad_state, double complex *state) {
    int i;
    non_rad_state[0*dim + 0] = gamma_nr * state[dim*1 + 1];
    for(i = 1; i < 3; i++) {
        non_rad_state[dim*i + i] = gamma_nr * (state[dim*(i+1) + i+1] - state[dim*i + i]);
    }
    non_rad_state[dim*3 + 3] = -gamma_nr * state[dim*3 + 3];
    non_rad_state[dim*4 + 4] = gamma_nr * state[dim*5 + 5];
    non_rad_state[dim*5 + 5] = gamma_nr * (state[dim*6 + 6] - state[dim*5 + 5]);
    non_rad_state[dim*6 + 6] = -gamma_nr * state[dim*6 + 6];
}

void SpontEm(double complex *spont_em_state, double complex *state) {
    int i;
    for(i = 0; i < 4; i++) {
        spont_em_state[dim*i + i] = spont_em * (state[dim*4 + 4] + state[dim*5 + 5] + state[dim*6 + 6]);
    }
    for(i = 4; i < dim; i++) {
        spont_em_state[dim*i + i] = -4 * spont_em * state[dim*i + i];
    }
}

/*void RK2_step(double complex *state, double dt, double temp, double j_0_3, double j_0_4) {
    int i,j;
    double complex *aux_state, *diff_state, *K1, *K2;
    aux_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    diff_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    K1 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K2 = (double complex*) calloc(dim*dim, sizeof(double complex));

    Diff(diff_state, state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            K1[dim*i + j] = diff_state[dim*i + j] * dt;
            aux_state[dim*i + j] = state[dim*i + j] + 0.5* K1[dim*i + j];
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            K2[dim*i + j] = diff_state[dim*i + j] * dt;
        }
    }
    for(i = 0; i < dim*dim; i++) {
        state[i] = state[i] + 0.5 * (K1[i] + K2[i]);
    }
    free(K1);
    free(K2);
    free(aux_state);
    free(diff_state);
}*/
/*
void RK4_step(double complex *state, double dt, double temp, double j_0_3, double j_0_4) {
    int i, j;
    
    double complex *aux_state, *diff_state, *K1, *K2, *K3, *K4;
    aux_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    diff_state = (double complex*) calloc(dim*dim, sizeof(double complex));
    K1 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K2 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K3 = (double complex*) calloc(dim*dim, sizeof(double complex));
    K4 = (double complex*) calloc(dim*dim, sizeof(double complex));

    Diff(diff_state, state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            K1[dim*i + j] = diff_state[dim*i + j] * dt;
            aux_state[dim*i + j] = state[dim*i + j] + 0.5* K1[dim*i + j];
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            K2[dim*i + j] = diff_state[dim*i + j] * dt;
            aux_state[dim*i + j] = state[dim*i + j] + 0.5* K2[dim*i + j];
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            K3[dim*i + j] = diff_state[dim*i + j] * dt;
            aux_state[dim*i + j] = state[dim*i + j] + K3[dim*i + j];
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            K4[dim*i + j] = diff_state[dim*i + j] * dt;
        }
    }
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            state[dim*i + j] += (double)1/6 * (K1[dim*i + j] + 2 * K2[dim*i + j] + 2 * K3[dim*i + j] + K4[dim*i + j]); 
        }
    }
    free(aux_state);
    free(diff_state);
    free(K1);
    free(K2);
    free(K3);
    free(K4);
} */