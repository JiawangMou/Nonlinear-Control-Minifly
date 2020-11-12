/*******************************************************************************
* @Copylift (c) 2020, Jiawang Mou, Inc.
* @
* @ pluginTemplate : [description]
* @
* @ filename : matrix.h
* @ author   : Jiawang Mou(moujiawang@sjtu.edu.cn)
* @ create   : 2020/11/11 	 10:45:43
******************************************************************************/

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
void MulMatrixDD(float *x,float *y, int m,int k,int n, float *z);
void TransSquareD(float *x, int m);
void JacbiCor(double * pMatrix,int nDim, double *pdblVects, double *pdbEigenValues, double dbEps,int nJt);
void TransMatrixD(float *x, int m, int n, float *z);
void TransMatrixS(short *x, int m, int n, float *z);
void MatrixADD(float *x,float *y,int m,int k, float *z);
void vector3_crossproduct(float *x,float *y, float *output);
#endif	 //__MATRIX_H__
////////////////////////////////// EOF /////////////////////////////////////////
