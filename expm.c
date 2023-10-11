/*
This algotihtm is based on the paper:
@article{ibanez2022two,
  title={Two taylor algorithms for computing the action of the matrix
exponential on a vector}, author={Ib{\'a}{\~n}ez, Javier and Alonso, Jos{\'e} M
and Alonso-Jord{\'a}, Pedro and Defez, Emilio and Sastre, Jorge},
  journal={Algorithms},
  volume={15},
  number={2},
  pages={48},
  year={2022},
  publisher={MDPI}
}
The second algorithm is implemented in this code.
*/

#include "functions.h"
#include "phg.h"
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <phg/map.h>
#include <phg/matvec.h>
#include <phg/phg-mpi.h>
#include <phg/utils.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#define check phgPrintf("HELLO I AM HERE LINE=%d\n", __LINE__)
#define Pi 3.14159265359
#define h_bar 1
#define mass 1
static int min(int a, int b) {
  if (a < b) {
    return a;
  } else {
    return b;
  }
}
//==================================compute
// Av=expm(A*deltat)*v==========================
void matrix_expontenial(MAT *A, VEC *v, FLOAT deltat, int m, int mmax, int smax,
                        double tol, VEC **Av) {

  MAT *Adeltat;
  Adeltat = phgMatAXPBY(deltat, A, 0.0, NULL);
  VEC *v1;
  v1 = phgMatVec(0.0, 1.0, Adeltat, v, 0.0, NULL);
  int j;
  VEC *v1_tmp = phgVecCopy(v1, NULL);
  for (j = 2; j < (m + 2); j++) {
    phgMatVec(0.0, 1.0, Adeltat, v1, 0.0, &v1_tmp);
    phgVecCopy(v1_tmp, &v1);
  }
  int s;
  double nfac[mmax];
  nfac[0] = 1;
  for (j = 1; j < mmax; j++) {

    nfac[j] = nfac[j - 1] / (j + 1);
  }
  s = (int)ceil(pow(phgVecNorm1(v1, 0, NULL) * nfac[m] / tol, 1. / (m + 1)));
  FLOAT v1norm;
  float powtmp = pow(phgVecNorm1(v1, 0, NULL) * nfac[m] / tol, 1. / (m + 1));
  double p = m * s;
  double s1, p1;

  int f = 0;
  while (f == 0 && m < mmax) {
    m++;
    phgMatVec(0.0, 1.0, Adeltat, v1, 0.0, &v1_tmp);
    phgVecCopy(v1_tmp, &v1);
    s1 = (int)ceil(pow(phgVecNorm1(v1, 0, NULL) * nfac[m] / tol, 1. / (m + 1)));

    p1 = m * s1;
    if (p1 <= p) {
      p = p1;
      s = s1;
    } else {
      m = m - 1;
      f = 1;
    }
  }
  s = min(s, smax);
  p = m * s;
  VEC *w;
  w = phgVecCopy(v, NULL);
  phgVecCopy(v, &v1);
  for (j = 0; j < m; j++) {
    phgMatVec(0, 1.0, Adeltat, v1, 0.0, &v1_tmp);
    phgVecCopy(v1_tmp, &v1);
    phgVecAXPBY(nfac[j] / pow(s, (j + 1)), v1, 1., &w);
  }

  MAT *Adeltats;
  Adeltats = phgMatAXPBY(1. / s, Adeltat, 0.0, NULL);
  int i;
  VEC *v2, *v2tmp;
  v2tmp = phgVecCopy(w, NULL);
  v2 = phgVecCopy(w, NULL);
  for (i = 0; i < (s - 1); i++) {

    phgVecCopy(w, &v2);
    for (j = 0; j < m; j++) {
      phgMatVec(0.0, 1.0, Adeltats, v2, 0.0, &v2tmp);
      phgVecCopy(v2tmp, &v2);
      phgVecAXPBY(nfac[j], v2, 1.0, &w);
    }
  }

  phgVecCopy(w, Av);
  phgMatDestroy(&Adeltats);
  phgVecDestroy(&v2);
  phgMatDestroy(&Adeltat);
  phgVecDestroy(&v1);
  phgVecDestroy(&w);
}
