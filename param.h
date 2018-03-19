#ifndef PARAM_H
#define PARAM_H

#define ROUND(x) (x<0?ceil((x)-0.5):floor((x)+0.5))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#define SICORTEX 0
#define DEBUG 1
#define DSTEP 1
#define PRINTID 0

#define PI 3.14159265

#define NX 150
#define NY NX
#define N_WOUND 300
#define MAXT 10

#define LX 1
#define LY 1

#define ANSYS_SCALING 100

#define A_ELL 0.25
#define B_ELL 0.25
#define PHI_ELL 0
#define WOUND_OPTION_ELL 1

#define P_MITOSIS 0.05
#define P_MITOSIS_INSIDE 0.02
#define Q_MOBILITY 0.1
#define Q_DIFFUSION 0.08
#define COEF_E0 0.6
#define COEF_E1 1.5
#define RATE_OF_DECAY_OF_GF 0

#define MAX_DIFFUSION_MOBILITY ROUND(Q_MOBILITY*MAX(NX, NY)) // here is the
#define MAX_DIFFUSION_GROWTH_FACTOR ROUND(Q_DIFFUSION*MAX(NX, NY)) // problem for scalability

double X_ansys_E[N_WOUND-2];
double Y_ansys_E[N_WOUND-2];

double Energy[N_WOUND-2];
double Energy_0;

//debug_debut
/*int index_rand;
int index_rand_tot;
double tab_rand[10000];
double rand_debug(void);*/
//debug_fin

#endif
