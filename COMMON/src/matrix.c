#include "math.h"
#include "matrix.h"

//矩阵乘法
/********参数表*******
@Parameter    x:    m行k列矩阵（用一维数组表示）
@Parameter    y:    k行n列矩阵（用一维数组表示）
@Parameter    m,k,n:    矩阵行列参数
@Parameter    z:    m行n列输出矩阵（用一维数组表示）
***********************/
void MulMatrixDD(float *x,float *y, int m,int k,int n, float *z)
{
    for(int nm=0; nm<m; nm++)
        for(int nn=0; nn<n; nn++)
            for(int nk=0; nk<k; nk++)
                z[nm*n+nn] += x[nm*k+nk]*y[nk*n+nn];

}

//矩阵加法
/********参数表*******
@Parameter    x:    m行k列矩阵（用一维数组表示）
@Parameter    y:    m行k列矩阵（用一维数组表示）
@Parameter    m,k:    矩阵行列参数
@Parameter    z:    m行k列输出矩阵（用一维数组表示）
***********************/
void MatrixADD(float *x,float *y,int m,int k, float *z)
{
    for(int i=0; i<m; i++)
        for(int j=0; j<k; j++)
            z[i*m+j] = x[i*m+j] + y[i*m+j];
}

//方阵转置
/********参数表*******
@Parameter    x:    m行m列矩阵（用一维数组表示）
@Parameter    m:    矩阵行列数
***********************/
void TransSquareD(float *x, int m)
{
    float temp;
    for(int nm=0; nm<m; nm++){            //对原矩阵第nm行
        for(int nn=0; nn<nm; nn++){        //对原矩阵第nn列
            temp = x[nm*m+nn];            //z矩阵第nn行第nm列
            x[nm*m+nn] = x[nn*m+nm];
            x[nn*m+nm] = temp;}}
}

//非方阵转置
/********参数表*******
@Parameter    x:    m行n列矩阵（用一维数组表示）
@Parameter    m,n:    矩阵行列数
@Parameter    z:    n行m列矩阵（用一维数组表示）
***********************/
void TransMatrixD(float *x, int m, int n, float *z)
{
    for(int nm=0; nm<m; nm++)            //对原矩阵第nm行
        for(int nn=0; nn<n; nn++)        //对原矩阵第nn列
            z[nn*m+nm] = x[nm*n+nn];    //z矩阵第nn行第nm列
}
void TransMatrixS(short *x, int m, int n, float *z)
{
    for(int nm=0; nm<m; nm++)            //对原矩阵第nm行
        for(int nn=0; nn<n; nn++)        //对原矩阵第nn列
            z[nn*m+nm] = (float)x[nm*n+nn];    //z矩阵第nn行第nm列
}

/********参数表*******
* @brief 求实对称矩阵的特征值及特征向量的雅克比法  
* 利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量  
@Parameter    pMatrix                长度为n*n的数组，存放实对称矩阵
@Parameter    nDim                   矩阵的阶数  
@Parameter    pdblVects              长度为n*n的数组，返回特征向量(按列存储)  
@Parameter    dbEps                  精度要求  
@Parameter    nJt                    整型变量，控制最大迭代次数  
@Parameter    pdbEigenValues         特征值数组 
***********************/
void JacbiCor(double * pMatrix,int nDim, double *pdblVects, double *pdbEigenValues, double dbEps,int nJt)  
{  
    int i,j;


    for(i = 0; i < nDim; i ++)   
    {     
        pdblVects[i*nDim+i] = 1.0f;   
        for(int j = 0; j < nDim; j ++)   
        {   
            if(i != j)     
                pdblVects[i*nDim+j]=0.0f;   
        }   
    }   
  
    int nCount = 0;     //迭代次数  
    while(1)  
    {  
        //在pMatrix的非对角线上找到最大元素  
        double dbMax = pMatrix[1];  
        int nRow = 0;  
        int nCol = 1;  
        for (i = 0; i < nDim; i ++)          //行  
        {  
            for (j = 0; j < nDim; j ++)      //列  
            {  
                double d = fabs(pMatrix[i*nDim+j]);   
  
                if((i!=j) && (d> dbMax))   
                {   
                    dbMax = d;     
                    nRow = i;     
                    nCol = j;   
                }   
            }  
        }  
  
        if(dbMax < dbEps)     //精度符合要求   
            break;    
  
        if(nCount > nJt)       //迭代次数超过限制  
            break;  
  
        nCount++;  
  
        double dbApp = pMatrix[nRow*nDim+nRow];  
        double dbApq = pMatrix[nRow*nDim+nCol];  
        double dbAqq = pMatrix[nCol*nDim+nCol];  
  
        //计算旋转角度  
        double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);  
        double dbSinTheta = sin(dbAngle);  
        double dbCosTheta = cos(dbAngle);  
        double dbSin2Theta = sin(2*dbAngle);  
        double dbCos2Theta = cos(2*dbAngle);  
  
        pMatrix[nRow*nDim+nRow] = dbApp*dbCosTheta*dbCosTheta +   
            dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;  
        pMatrix[nCol*nDim+nCol] = dbApp*dbSinTheta*dbSinTheta +   
            dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;  
        pMatrix[nRow*nDim+nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;  
        pMatrix[nCol*nDim+nRow] = pMatrix[nRow*nDim+nCol];  
  
        for(i = 0; i < nDim; i ++)   
        {   
            if((i!=nCol) && (i!=nRow))   
            {   
                int u = i*nDim + nRow;  //p    
                int w = i*nDim + nCol;  //q  
                dbMax = pMatrix[u];   
                pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta;   
                pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta;   
            }   
        }   
  
        for (j = 0; j < nDim; j ++)  
        {  
            if((j!=nCol) && (j!=nRow))   
            {   
                int u = nRow*nDim + j;  //p  
                int w = nCol*nDim + j;  //q  
                dbMax = pMatrix[u];   
                pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta;   
                pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta;   
            }   
        }  
  
        //计算特征向量  
        for(i = 0; i < nDim; i ++)   
        {   
            int u = i*nDim + nRow;      //p     
            int w = i*nDim + nCol;      //q  
            dbMax = pdblVects[u];   
            pdblVects[u] = pdblVects[w]*dbSinTheta + dbMax*dbCosTheta;   
            pdblVects[w] = pdblVects[w]*dbCosTheta - dbMax*dbSinTheta;   
        }   
  
    }  
    for(i = 0; i < nDim; i ++)   
    {     
        pdbEigenValues[i] = pMatrix[i*nDim+i];  
    }   

    //设定正负号  
    for(i = 0; i < nDim; i ++)   
    {  
        double dSumVec = 0;  
        for(j = 0; j < nDim; j ++)  
            dSumVec += pdblVects[j * nDim + i];  
        if(dSumVec<0)  
        {  
            for(j = 0;j < nDim; j ++)  
                pdblVects[j * nDim + i] *= -1;  
        }  
    } 

}


/********************************向量相关运算*******************************************/
void vector3_crossproduct(float *x,float *y, float *output)
{
    output[0] = x[1] * y[2] - x[2] * y[1];
    output[1] = x[2] * y[0] - x[0] * y[2];
    output[2] = x[0] * y[1] - x[1] * y[0];
}
