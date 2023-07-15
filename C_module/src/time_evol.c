#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "aux.h"
#include "time_evol.h"
void Time_evol(double complex *state, double tf, double dt, double temp, double j_0_3, double j_0_4){
    int N = (int) (tf/dt), i;
    printf("N = %d", N);
    for(i = 0; i < N; i++){
        RK4_step(state, dt, temp, j_0_3, j_0_4);
        printf("\r Progress = %.2lf", (double)i/N*100);
    }
}

void RK4_step(double complex *state, double dt, double temp, double j_0_3, double j_0_4){
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
            *(K1 + dim*i + j) = *(diff_state + dim*i + j) * dt;
            *(aux_state + dim*i + j) = *(state + dim*i + j) + 0.5* *(K1 + dim*i + j);
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            *(K2 + dim*i + j) = *(diff_state + dim*i + j) * dt;
            *(aux_state + dim*i + j) = *(state + dim*i + j) + 0.5* *(K2 + dim*i + j);
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            *(K3 + dim*i + j) = *(diff_state + dim*i + j) * dt;
            *(aux_state + dim*i + j) = *(state + dim*i + j) + *(K3 + dim*i + j);
        }
    }
    Diff(diff_state, aux_state, temp, j_0_3, j_0_4);    
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            *(K4 + dim*i + j) = *(diff_state + dim*i + j) * dt;
        }
    }
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            *(state + dim*i + j) += (double)1/6 * (*(K1 + dim*i + j) + 2 * *(K2 + dim*i + j) + 2 * *(K3 + dim*i + j) * *(K4 + dim*i + j)); 
        }
    }
    free(aux_state);
    free(diff_state);
    free(K1);
    free(K2);
    free(K3);
    free(K4);
}

void Diff(double complex *diff_state, double complex *state, double temp, double j_0_3, double j_0_4){
    int i, j;
    double Rabi4;
    double *ks, *ws;
    ks = (double *) calloc(6, sizeof(double));
    ws = (double *) calloc(6, sizeof(double));
    Rabi4 = Rabifreq(j_0_4);
    *ks = k12; *(ks + 1) = k23; *(ks + 2) = k34; *(ks + 3) = (double)Rabi4/sqrt(Thermal_num(temp, w)); *(ks + 4) = k56; *(ks + 5) = k67;  
    *ws = w1; *(ws + 1) = w2; *(ws + 2) = w3; *(ws + 3) = w; *(ws + 4) = w5; *(ws + 5) = w6; 

    for(i = 1; i < dim-1; i++){
        for(j = 1; j < dim-1; j++){
            *(diff_state + dim*i + j) = 1*I * *(state + dim*i + j + 1) * sqrt(Thermal_num(temp, *(ws + j))) * *(ks+j) 
            + 1*I * *(state + dim*i + j - 1) * sqrt(Thermal_num(temp, *(ws + j-1))) * *(ks+j-1)
            - 1*I * *(state + dim*(i+1) + j) * sqrt(Thermal_num(temp, *(ws + i))) * *(ks + i)
            - 1*I * *(state + dim*(i-1) + j) * sqrt(Thermal_num(temp, *(ws + i-1))) * *(ks + i-1);
        }
    }
    for(j = 1; j < dim-1; j++){
        *(diff_state + dim*0 + j) = 1*I * *(state + dim*0 + j + 1) * sqrt(Thermal_num(temp, *(ws + j))) * *(ks+j) 
        + 1*I * *(state + dim*0 + j - 1) * sqrt(Thermal_num(temp, *(ws + j-1))) * *(ks+j-1)
        - 1*I * *(state + dim*(0+1) + j) * sqrt(Thermal_num(temp, *(ws + 0))) * *(ks + 0);
    }
    for(j = 1; j < dim-1; j++){
        *(diff_state + dim*(dim-1) + j) = 1*I * *(state + dim*(dim-1) + j + 1) * sqrt(Thermal_num(temp, *(ws + j))) * *(ks+j) 
        + 1*I * *(state + dim*(dim-1) + j - 1) * sqrt(Thermal_num(temp, *(ws + j-1))) * *(ks+j-1)
        - 1*I * *(state + dim*((dim-1)-1) + j) * sqrt(Thermal_num(temp, *(ws + (dim-1)-1))) * *(ks + (dim-1)-1);
    }
    for(i = 1; i < dim-1; i++){
        *(diff_state + dim*i + 0) = 1*I * *(state + dim*i + 0 + 1) * sqrt(Thermal_num(temp, *(ws + 0))) * *(ks+0) 
        - 1*I * *(state + dim*(i+1) + 0) * sqrt(Thermal_num(temp, *(ws + i))) * *(ks + i)
        - 1*I * *(state + dim*(i-1) + 0) * sqrt(Thermal_num(temp, *(ws + i-1))) * *(ks + i-1);
    }
    for(i = 1; i < dim-1; i++){
        *(diff_state + dim*i + (dim-1)) = 1*I * *(state + dim*i + (dim-1) - 1) * sqrt(Thermal_num(temp, *(ws + (dim-1)-1))) * *(ks+(dim-1)-1)
        - 1*I * *(state + dim*(i+1) + (dim-1)) * sqrt(Thermal_num(temp, *(ws + i))) * *(ks + i)
        - 1*I * *(state + dim*(i-1) + (dim-1)) * sqrt(Thermal_num(temp, *(ws + i-1))) * *(ks + i-1);
    } 
    *(diff_state + dim*0 + 0) = 1*I * *(state + dim*0 + 0 + 1) * sqrt(Thermal_num(temp, *(ws + 0))) * *(ks+0) 
    - 1*I * *(state + dim*(0+1) + 0) * sqrt(Thermal_num(temp, *(ws + 0))) * *(ks + 0);

    *(diff_state + dim*(dim-1) + (dim-1)) = 1*I * *(state + dim*(dim-1) + (dim-1) - 1) * sqrt(Thermal_num(temp, *(ws + (dim-1)-1))) * *(ks+(dim-1)-1)
    - 1*I * *(state + dim*((dim-1)-1) + (dim-1)) * sqrt(Thermal_num(temp, *(ws + (dim-1)-1))) * *(ks + (dim-1)-1);

    for(i = 0; i < dim-1; i++){
        if(i != 3){
            *(diff_state + dim*i + i) += gamma_nr * *(state + dim*(i+1) + (i+1));
        }
    }    
    for(i = 1; i < dim; i++){
        if(i != 4){
            *(diff_state + dim*i + i) -= gamma_nr * *(state + dim*i + i);
        }
    }

    for(i = 0; i < 4; i++){
        *(diff_state + dim*i + i) += spont_em*(*(state + dim*4 + 4) + *(state + dim*5 + 5) + *(state + dim*6 + 6));
    }    
    for(i = 4; i < dim; i++){
        *(diff_state + dim*i + i) -= 4 * spont_em * *(state + dim*i + i);
    }

    for(i = 0; i < 3; i++){
        for(j = 4; j < dim; j++){
            *(diff_state + dim*i + j) -= 0.5 * gamma_nr * *(state + dim*i + j);
        }
    }
    for(i = 4; i < dim; i++){
        for(j = 0; j < 3; j++){
            *(diff_state + dim*i + j) -= 0.5 * gamma_nr * *(state + dim*i + j);
        }
    }

    for(i = 0; i < dim-1; i++){
        *(diff_state + dim*i + (i+1)) -= 0.5*gamma_nr * *(state + dim*i + i+1);
        *(diff_state + dim*(i+1) + i) -= 0.5*gamma_nr * *(state + dim*(i+1) + i);
    }
    for(i = 1; i< dim; i++){
        *(diff_state + dim*i + (i-1)) -= 0.5*gamma_nr * *(state + dim*i + i-1);
        *(diff_state + dim*(i-1) + i) -= 0.5*gamma_nr * *(state + dim*(i-1) + i);
    }


    free(ks);
    free(ws);
}
