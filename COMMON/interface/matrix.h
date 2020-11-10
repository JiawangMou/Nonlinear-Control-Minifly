#ifndef __MATRIX_H
#define __MATRIX_H

void MulMatrixDD(double *x,double *y, int m,int k,int n, double *z);
void TransSquareD(double *x, int m);
void JacbiCor(double * pMatrix,int nDim, double *pdblVects, double *pdbEigenValues, double dbEps,int nJt);
void TransMatrixD(double *x, int m, int n, double *z);
void TransMatrixS(double *x, int m, int n, double *z);
#endif