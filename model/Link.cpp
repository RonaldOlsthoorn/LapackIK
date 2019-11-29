//
// Created by ronald on 26-11-19.
//

#include <cstdio>
#include <iostream>
#include "Link.h"
#include "math.h"
#include "mkl.h"

double *Link::getForwardMatrix() {

    printf("radius: %f \n", radius);

    double *A = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    std::cout<< "A adress: "<<A;
    A[0] = cos(jointAngle);
    A[1] = -sin(jointAngle);
    A[4] = sin(jointAngle);
    A[5] = cos(jointAngle);
    A[10] = 1;
    A[15] = 1;

    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            printf ("%12.5G", A[4*i+j]);
        }
        printf("\n");
    }



    double *B = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    B[0] = 1;
    B[5] = 1;
    B[10] = 1;
    B[11] = offset;
    B[15] = 1;

    double *C = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    C[0] = 1;
    C[5] = 1;
    C[10] = 1;
    C[3] = radius;
    C[15] = 1;

    double *D = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    D[0] = 1;
    D[5] = cos(twist);
    D[6] = -sin(twist);
    D[9] = sin(twist);
    D[10] = cos(twist);
    D[15] = 1;

    double *scrap1 = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 4, 4, 1, A, 4, B, 4, 0, scrap1, 4);

    double *scrap2 = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 4, 4, 1, scrap1, 4, C, 4, 0, scrap2, 4);

    double *res = (double *)mkl_malloc( 4*4*sizeof( double ), 64);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 4, 4, 1, scrap2, 4, D, 4, 0, res, 4);



    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            printf ("%12.5G", B[4*i+j]);
        }
        printf("\n");
    }

    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            printf ("%12.5G", C[4*i+j]);
        }
        printf("\n");
    }

    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            printf ("%12.5G", D[4*i+j]);
        }
        printf("\n");
    }

    mkl_free(A);
    mkl_free(B);
    mkl_free(C);
    mkl_free(D);
    mkl_free(scrap1);
    mkl_free(scrap2);

    return res;
}
