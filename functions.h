#include "phg.h"
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <string.h>
// define the domain of inclusion in a unit cell: [L1, L2]^{3}
#define L1 0.25
#define L2 0.75

// Number of inclusion in each direction
#define noc 8

// the relative size of periodic microstructure
#define Epsilon 0.125

// the number of electrons
#define NEV 4

// the parameters used in the computation of electron density and current density
#define E_f 80
#define KT 1
#define Nd 0

// parameters controling the nonlinear iteration
#define alpha 1.0
#define rtol 5.0e-9

// if using the exchange-correlation potential, then set xc=1
#define xc 0

// the value of permeability 
#define  in_muz 1     //inclusion
#define  out_muz 200  //matrix

// the value of effective mass
#define  in_indexm 0.005   //inclusion
#define  out_indexm 2  //matrix

// the value of confining potential
#define  in_indexv  4    //inclusion
#define  out_indexv 1  //matrix


#define nsteps 501//201//1001
#define stime 0.0025
 

// number of basis
//#define Nbasis 42925
//91
//number of Energy
//#define Energynumber 100//6



extern FLOAT ctime;

extern FLOAT ptime;
extern int   global_flag;

//extern void cell_maxwell(DOF *w_h1, DOF *w_h2, DOF *w_h3, DOF *n_h1, DOF *n_h2, DOF *n_h3, double **mat_coef);

//extern void maxwell_addcorrector(DOF *u_h, DOF *v_h2, DOF *v_hf, DOF *w_h11, DOF *w_h21, DOF *w_h31, DOF *n_h11, DOF *n_h21, DOF *n_h31);

//extern void  maxwell_solver(DOF *u_h, DOF *u_p, DOF *u_pp, DOF *H);



//extern void solve_eigen_problem(DOF **u_h, FLOAT *evals, DOF *emassmat, DOF *permittivity, DOF *vconf);

//extern FLOAT cell_schrodinger(DOF * N1_h1, DOF * N2_h1, DOF * N3_h1, DOF * N11_h1,
  //                DOF * N12_h1, DOF * N13_h1, DOF * N21_h1, DOF * N22_h1,
    //              DOF * N23_h1, DOF * N31_h1, DOF * N32_h1, DOF * N33_h1, FLOAT * v0);

//extern void schrod_addcorrector(DOF *u0_h, DOF *u1_h, DOF *u2_h, DOF *Nu1_h, DOF *Nu2_h,
  //             DOF *Nu3_h, DOF *Nu11_h, DOF *Nu12_h, DOF *Nu13_h, DOF *Nu21_h,
    //           DOF *Nu22_h, DOF *Nu23_h, DOF *Nu31_h, DOF *Nu32_h, DOF *Nu33_h);

//extern void schrodinger_solver(DOF *u_re, DOF *u_im, DOF *u_pre, DOF *u_pim, DOF *efield, DOF *emassmat, DOF *vconf, DOF *vxc);

//extern void calculate_u1minusu0(DOF *u1_h, DOF *u0_h, DOF *N1_h, DOF *N2_h, DOF *N3_h);

#define bzero(buffer, size) memset(buffer, 0, size)
/*
extern void CN_solver_parallel(FLOAT t_ini, VEC **c, VEC *cini,
                        INT tN, FLOAT deltat, FLOAT mass, FLOAT *E,
                        MAT *DiagMat, MAT *x1Mat_imag, MAT *x2Mat_imag,
                        MAT *x3Mat_imag, MAT *I2n);*/
extern void CN_solver_parallel(FLOAT t_ini, VEC **c, VEC *cini, INT tN, FLOAT deltat,
                        FLOAT mass, MAT *Atilde, MAT *I2n);

extern void matrix_expontenial(MAT *A, VEC *v, FLOAT deltat, int m, int mmax, int smax,
                        double tol, VEC **Av);
