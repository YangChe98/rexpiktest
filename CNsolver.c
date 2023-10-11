
#include "phg.h"
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <phg/map.h>
#include <phg/matvec.h>
#include <phg/phg-mpi.h>
#include <phg/solver.h>
#include <phg/utils.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#define check phgPrintf("HELLO I AM HERE LINE=%d\n", __LINE__)
#define Pi 3.14159265359
#define h_bar 1
/*
void CN_solver_parallel(FLOAT t_ini, VEC **c, VEC *cini, INT tN, FLOAT deltat,
                        FLOAT mass, FLOAT *E, MAT *DiagMat, MAT *x1Mat_imag,
                        MAT *x2Mat_imag, MAT *x3Mat_imag, MAT *I2n) {*/

void CN_solver_parallel(FLOAT t_ini, VEC **c, VEC *cini, INT tN, FLOAT deltat,
                        FLOAT mass,  MAT *Atilde, MAT *I2n) {
  VEC *c_pre;

  FLOAT t_now;
  FLOAT i;
  t_now = t_ini;

//  check;
  c_pre = phgVecCopy(cini, NULL);
//  check;
/*
 MAT *Atilde;
//  check;
  Atilde = phgMatAXPBY(1.0, DiagMat, 0, NULL);
//  check;
  phgMatAXPBY(E[0], x1Mat_imag, 1.0, &Atilde);
  phgMatAXPBY(E[1], x2Mat_imag, 1.0, &Atilde);
  phgMatAXPBY(E[2], x3Mat_imag, 1.0, &Atilde);
  */
//  check;
  MAT *Ahat, *Ahat2;
  Ahat = phgMatAXPBY(h_bar, I2n, 0, NULL);
  Ahat2 = phgMatAXPBY(h_bar, I2n, 0, NULL);
  phgMatAXPBY(-deltat / 2, Atilde, 1.0, &Ahat);
  phgMatAXPBY(deltat / 2, Atilde, 1.0, &Ahat2);
  //  phgMatDumpMATLAB(Ahat, "Ahat1", "Ahat11.m");
  // check;
  SOLVER *matsolver;
// matsolver=phgSolverMat2Solver(SOLVER_SUPERLU,Ahat); 
 //matsolver = phgSolverMat2Solver(SOLVER_PETSC, Ahat);
  matsolver = phgSolverMat2Solver(SOLVER_DEFAULT, Ahat);
  phgSolverSetMaxIt(matsolver, 1e6);
  phgSolverSetTol(matsolver, 1e-12);
  // check;

  VEC *rhss;

  MAT *solvermat;
  solvermat = phgSolverGetMat(matsolver);
  //  phgMatDumpMATLAB(solvermat, "solvermat1", "solvermat.m");
 // check;
  VEC *x;
  phgSolverAssemble(matsolver);
  for (i = 0; i < tN; i++) {
    phgSolverResetRHS(matsolver);
    rhss = phgMatVec(0, 1.0, Ahat2, c_pre, 0.0, NULL);
    x = phgMapCreateVec(matsolver->mat->rmap, 1);
    matsolver->rhs->data = rhss->data;
    matsolver->rhs->assembled = TRUE;
    matsolver->rhs_updated = TRUE;

    // phgVecDumpMATLAB(matsolver->rhs, "rhss1", "rhss.m");
    // check;
    phgSolverVecSolve(matsolver, FALSE, x);
    //    if (i == 1) {
    //    phgVecDumpMATLAB(matsolver->rhs, "rhs11", "rhs1.m");
    //}
  //  phgPrintf("schrodinger nits=%d, resid=%0.4lg \n", matsolver->nits,
  //            (double)matsolver->residual);

   // phgPrintf("schrodinger norm2= %e \n", phgVecNorm2(x, 0, NULL));
   // phgPrintf("schrodinger norm2= %e \n", phgVecNorm2(x, 0, NULL));
    //    phgVecDumpMATLAB(c_pre, "c_pre1", "c_pre.m");
    phgVecCopy(x, &c_pre);
    // if (i == 0) {
    //    phgVecDumpMATLAB(c_pre, "c_pre2", "c_pre3.m");
    //}

    //    phgVecDumpMATLAB(x, "x", "x1.m");
  }
  phgVecCopy(x, c);
  //phgMatDestroy(&Atilde);
  phgMatDestroy(&Ahat);
  phgMatDestroy(&Ahat2);
  phgMatDestroy(&solvermat);
  phgSolverDestroy(&matsolver);
  // phgVecDumpMATLAB(c, "c1", "c.m");
}
