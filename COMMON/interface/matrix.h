/*******************************************************************************
* @Copylift (c) 2020, Jiawang Mou, Inc.
* @
* @ pluginTemplate : [description]
* @
* @ filename : matrix.h
* @ author   : Jiawang Mou(moujiawang@sjtu.edu.cn)
* @ create   : 2020/11/11 	 10:45:43
/******************************************************************************/

#ifndef __MATRIX_H__
#define __MATRIX_H__

////////////////////////////////////////////////////////////////////////////////
// Headers
//

////////////////////////////////////////////////////////////////////////////////
// Typedefs & Constants
//

////////////////////////////////////////////////////////////////////////////////
// Classes
//

////////////////////////////////////////////////////////////////////////////////
// Functions
//
void MulMatrixDD(double *x,double *y, int m,int k,int n, double *z);
void TransSquareD(double *x, int m);
void JacbiCor(double * pMatrix,int nDim, double *pdblVects, double *pdbEigenValues, double dbEps,int nJt);
void TransMatrixD(double *x, int m, int n, double *z);
void TransMatrixS(double *x, int m, int n, double *z);
void MatrixADD(float *x,float *y,int m,int k, float *z);
void vector3_crossproduct(float *x,float *y, float *output);
#endif	 //__MATRIX_H__
////////////////////////////////// EOF /////////////////////////////////////////