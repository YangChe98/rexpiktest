#include "phg.h"
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <phg/matvec.h>
#include <phg/phg-mpi.h>
#include <phg/utils.h>
#include <stdlib.h>

#define check phgPrintf("HELLO I AM HERE LINE=%d\n", __LINE__)
#define Pi 3.14159265359
#define h_bar 1
#define mass 1

int main(int argc, char *argv[]) {
  phgOptionsPreset("-solver_rtol 1.0e-9");
  phgInit(&argc, &argv);

  INT Nbasis = 40;
  INT Nbasis2 = Nbasis * Nbasis;
  INT Nbasis3 = Nbasis * Nbasis * Nbasis;

  INT nprocs = phgNProcs;
  INT vec1procs = 0;
  if (2 * Nbasis * Nbasis * Nbasis % nprocs == 0) {
    vec1procs = 2 * Nbasis * Nbasis * Nbasis / nprocs;
  } else {
    return 0;
  }
  MAP *vecmap =
      phgMapCreateSimpleMap(phgComm, vec1procs, 2 * Nbasis * Nbasis * Nbasis);

  int m = 50;
  int mmax = 100;
  int smax = 100;
  double tol = 1e-16;

  VEC *cini = phgMapCreateVec(vecmap, 1);

  // VEC *c1 = NULL;
  // VEC *c2 = NULL;
  // VEC *c0 = NULL;
  VEC *c = NULL;

  //  VEC *Ac1 = NULL;
  // VEC *Ac2 = NULL;
  // VEC *Ac0 = NULL;
  VEC *Ac = NULL;
  MAT *Atilde = NULL;
  //  MAT *Atilde0 = NULL;
  // MAT *Atilde1 = NULL;
  // MAT *Atilde2 = NULL;
  phgVecDisassemble(cini);

  INT i;

  /*======================================MATRIX DIOPLE
   * MATRIX=======================================================*/
  VEC *c_pre;  // *c1_pre, *c0_pre, *c2_pre;
  VEC *Ac_pre; //*Ac1_pre, *Ac2_pre, *Ac0_pre;
  FLOAT x3value, x2value, x1value;
  INT I[vec1procs], I1[vec1procs], I2[2 * Nbasis3];

  // check;
  INT j, l, lll;
  FLOAT Ivalue[vec1procs];
  // check;
  FLOAT Diagvalue, I2nvalue;
  // check;
  FLOAT x1valuemat[2 * Nbasis3];
  // check;
  FLOAT x2valuemat[2 * Nbasis3];
  // check;
  FLOAT x3valuemat[2 * Nbasis3];

  MAT *DiagMat = phgMapCreateMat(vecmap, vecmap);
  MAT *x1Mat_imag = phgMapCreateMat(vecmap, vecmap);
  MAT *x2Mat_imag = phgMapCreateMat(vecmap, vecmap);
  MAT *x3Mat_imag = phgMapCreateMat(vecmap, vecmap);
  MAT *I2n = phgMapCreateMat(vecmap, vecmap);
  phgMatDisassemble(DiagMat);
  phgMatDisassemble(x1Mat_imag);
  phgMatDisassemble(x2Mat_imag);
  phgMatDisassemble(x3Mat_imag);
  phgMatDisassemble(I2n);

  INT i1, j1, k1, jjj, i2, j2, k2, j3, loctmp1, loctmp2, h;
  FLOAT ii, jj, kk;

  phgPrintf("%d \n.", vec1procs);
  for (l = 0; l < 2 * Nbasis3; l++) {

    I2[l] = l;
  }

  for (i = 0; i < vec1procs; i++) {
    jjj = DiagMat->rmap->partition[DiagMat->rmap->rank] + i;
    I[i] = jjj;
    I1[i] = (jjj + Nbasis3) % (2 * Nbasis3);

    j = jjj % Nbasis3;
    loctmp1 = (INT)jjj / Nbasis3;
    i1 = (INT)j / Nbasis2;
    j1 = (INT)(j - i1 * Nbasis2) / Nbasis;
    k1 = (j - i1 * Nbasis2 - j1 * Nbasis);
    ii = i1 - Nbasis / 2.;
    jj = j1 - Nbasis / 2.;
    kk = k1 - Nbasis / 2.;

    x3value = (ii * ii + jj * jj + kk * kk) / (2 * mass);

    if (jjj < Nbasis3) {
      Diagvalue = x3value;
    } else {
      Diagvalue = -x3value;
    }
    I2nvalue = 1;

    Ivalue[i] = j + 1.;

    for (l = 0; l < 2 * Nbasis3; l++) {
      lll = l;
      j3 = lll % Nbasis3;
      loctmp2 = (INT)lll / Nbasis3;
      i2 = (INT)j3 / Nbasis2;
      j2 = (INT)(j3 - i2 * Nbasis2) / Nbasis;
      k2 = (j3 - i2 * Nbasis2 - j2 * Nbasis);

      if (loctmp1 == loctmp2) {
        if (i1 == i2) {
          if (j1 == j2) {
            h = k1 + k2 - Nbasis;
            if (h == 0) {
              x3value = 0;
            } else {
              x3value = 8 * Pow(Pi, 3.0) * Pow(-1.0, h + 1) / (h);
            }
            x3valuemat[l] = x3value;
          } else {
            x3valuemat[l] = 0;
          }
        } else {
          x3valuemat[l] = 0;
        }
        if (i1 == i2) {
          if (k1 == k2) {
            h = j1 + j2 - Nbasis;
            if (h == 0) {
              x3value = 0;
            } else {
              x3value = 8 * Pow(Pi, 3.0) * Pow(-1.0, h + 1) / (h);
            }
            x2valuemat[l] = x3value;
          } else {
            x2valuemat[l] = 0;
          }
        } else {
          x2valuemat[l] = 0;
        }
        if (j1 == j2) {
          if (k1 == k2) {
            h = i1 + i2 - Nbasis;
            if (h == 0) {
              x3value = 0;
            } else {
              x3value = 8 * Pow(Pi, 3.0) * Pow(-1.0, h + 1) / (h);
            }
            x1valuemat[l] = x3value;
          } else {
            x1valuemat[l] = 0;
          }
        } else {
          x1valuemat[l] = 0;
        }
      } else {
        x3valuemat[l] = 0;
        x2valuemat[l] = 0;
        x1valuemat[l] = 0;
      }
    }

    phgMatAddGlobalEntry(DiagMat, I[i], I1[i], Diagvalue);
    phgMatAddGlobalEntry(I2n, I[i], I[i], I2nvalue);
    phgMatAddGlobalEntries(x1Mat_imag, 1, I + i, 2 * Nbasis3, I2, x1valuemat);
    phgMatAddGlobalEntries(x2Mat_imag, 1, I + i, 2 * Nbasis3, I2, x2valuemat);
    phgMatAddGlobalEntries(x3Mat_imag, 1, I + i, 2 * Nbasis3, I2, x3valuemat);
  }

  phgVecAddGlobalEntries(cini, 0, vec1procs, I, Ivalue);
  FLOAT cininorm;
  phgVecAssemble(cini);
  phgVecNorm2(cini, 0, &cininorm);
  phgPrintf("%e \n", cininorm);
  cini = phgVecAXPBY(1 / cininorm, cini, 0., NULL);
  phgMatAssemble(DiagMat);
  phgMatAssemble(x1Mat_imag);
  phgMatAssemble(x2Mat_imag);
  phgMatAssemble(x3Mat_imag);
  phgMatAssemble(I2n);

  int k;
  FILE *r_file;
  r_file = fopen("r_real.txt", "r");
  float *r_real;
  r_real = (float *)calloc((Nbasis + 1) * (Nbasis + 1) * (Nbasis + 1),
                           sizeof(FLOAT));
  if (phgRank == 0) {
    for (j = 0; j < (Nbasis + 1) * (Nbasis + 1); j++) {
      for (k = 0; k < (Nbasis + 1); k++) {
        i1 = j / (Nbasis + 1);
        j1 = j % (Nbasis + 1);
        fscanf(r_file, "%e\t",
               (r_real + i1 * (Nbasis + 1) * (Nbasis + 1) + j1 * (Nbasis + 1) +
                k));
      }
      fscanf(r_file, "\n ");
    }
  }
  fclose(r_file);
  FILE *output_file;
  output_file = fopen("output.txt", "w+");
  if (phgRank == 0) {
    for (i = 0; i < Nbasis + 1; i++) {
      for (j = 0; j < Nbasis + 1; j++) {
        for (k = 0; k < Nbasis + 1; k++) {
          fprintf(output_file, "%e\t",
                  *((r_real + i * (Nbasis + 1) * (Nbasis + 1) +
                     j * (Nbasis + 1) + k)));
        }
        fprintf(output_file, "\n");
      }
    }
  }
  fclose(output_file);
  phgFinalize();
  return 0;
}
