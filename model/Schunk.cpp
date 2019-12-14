//
// Created by ronald on 14-12-19.
//

#include <mkl_service.h>
#include <mkl.h>
#include <cstdio>
#include "Schunk.h"

void Schunk::getForwardMatrix(double *out) {

    links->front()->getForwardMatrix(out);

    printf("First joint \n");

    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            printf ("%12.5G", out[4*i+j]);
        }
        printf("\n");
    }

    for(std::vector<Link*>::iterator it = ++(links->begin()); it != links->end(); it++){
        double* forwardMatrix = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );
        (*it)->getForwardMatrix(forwardMatrix);

        printf("joint matrix \n");

        for(int i = 0; i<4; i++){
            for(int j = 0; j<4; j++){
                printf ("%12.5G", forwardMatrix[4*i+j]);
            }
            printf("\n");
        }

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    4, 4, 4, 1, out, 4, forwardMatrix, 4, 0, out, 4);

        printf("running matrix \n");
        for(int i = 0; i<4; i++){
            for(int j = 0; j<4; j++){
                printf ("%12.5G", out[4*i+j]);
            }
            printf("\n");
        }

        mkl_free(forwardMatrix);
    }
}
