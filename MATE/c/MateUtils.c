#include<stdio.h>
#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "mex.h"      /* MEX-file interface mechanism */
#include "mat.h"
#endif
void printmat(double *A, int nb_row, int nb_col)
{
    printf("\n");
    for (int i=0;i<nb_row; i++){
        for (int j=0;j<nb_col; j++){
            printf("%16.14e  ", A[i*nb_col+j]);
        }
        printf("\n");
    }
    printf("\n");
}
void printmatMEX(double *A, int nb_row, int nb_col)
{
    mexPrintf("\n");
    for (int i=0;i<nb_row; i++){
        for (int j=0;j<nb_col; j++){
            mexPrintf("%16.14e  ", A[i*nb_col+j]);
        }
        mexPrintf("\n");
    }
    mexPrintf("\n");
}
void printmatMAT(double *A, int nb_row, int nb_col, char name[10])
{
    mexPrintf("%s=[ \n",name);
    for (int i=0;i<nb_row; i++){
        for (int j=0;j<nb_col; j++){
            mexPrintf("%16.14e  ", A[i*nb_col+j]);
        }
        mexPrintf(";\n");
    }
    mexPrintf("]\n");
}