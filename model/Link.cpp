//
// Created by ronald on 26-11-19.
//

#include <cstdio>
#include <iostream>
#include "Link.h"
#include "math.h"
#include "mkl.h"

void Link::getForwardMatrix(double * out) {

    double *A = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    A[0] = cos(jointAngle);
    A[1] = -sin(jointAngle);
    A[3] = 0;
    A[4] = sin(jointAngle);
    A[5] = cos(jointAngle);
    A[6] = 0;
    A[7] = 0;
    A[8] = 0;
    A[9] = 0;
    A[10] = 1;
    A[11] = 0;
    A[12] = 0;
    A[13] = 0;
    A[14] = 0;
    A[15] = 1;

    double *B = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    B[0] = 1;
    B[2] = 0;
    B[3] = 0;
    B[4] = 0;
    B[5] = 1;
    B[6] = 0;
    B[7] = 0;
    B[8] = 0;
    B[9] = 0;
    B[10] = 1;
    B[11] = offset;
    B[12] = 0;
    B[13] = 0;
    B[14] = 0;
    B[15] = 1;

    double *C = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    C[0] = 1;
    C[1] = 0;
    C[2] = 0;
    C[3] = radius;
    C[4] = 0;
    C[5] = 1;
    C[6] = 0;
    C[7] = 0;
    C[8] = 0;
    C[9] = 0;
    C[10] = 1;
    C[11] = 0;
    C[12] = 0;
    C[13] = 0;
    C[14] = 0;
    C[15] = 1;

    double *D = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    D[0] = 1;
    D[1] = 0;
    D[2] = 0;
    D[3] = 0;
    D[4] = 0;
    D[5] = cos(twist);
    D[6] = -sin(twist);
    D[7] = 0;
    D[8] = 0;
    D[9] = sin(twist);
    D[10] = cos(twist);
    D[11] = 0;
    D[12] = 0;
    D[13] = 0;
    D[14] = 0;
    D[15] = 1;

    double *scrap1 = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 4, 4, 1, A, 4, B, 4, 0, scrap1, 4);

    double *scrap2 = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 4, 4, 1, scrap1, 4, C, 4, 0, scrap2, 4);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 4, 4, 1, scrap2, 4, D, 4, 0, out, 4);

    mkl_free(A);
    mkl_free(B);
    mkl_free(C);
    mkl_free(D);
    mkl_free(scrap1);
    mkl_free(scrap2);
}
