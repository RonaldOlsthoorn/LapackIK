//
// Created by ronald on 19-12-19.
//

#include <mkl.h>

void getRotationMatrix(double* fromUnit, double* toUnit, double *out){

    double *v = (double *)mkl_malloc( 1*3*sizeof( double ), 64 );
    double *w = (double *)mkl_malloc( 1*3*sizeof( double ), 64 );

    double dot = cblas_ddot (3, fromUnit, 1, toUnit, 1);

    for(int i=0;i<3;i++){
        v[i] = toUnit[i] - dot * fromUnit[i];
    }

    double length = cblas_dnrm2(3, v, 1);

    for(int i=0;i<3;i++){
        v[i] = v[i]/length;
    }

    double *fromMatrix = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );
    fromMatrix[0] = 0;
    fromMatrix[1] = -fromUnit[2];
    fromMatrix[2] = fromUnit[1];
    fromMatrix[3] = fromUnit[2];
    fromMatrix[4] = 0;
    fromMatrix[5] = -fromUnit[0];
    fromMatrix[6] = -fromUnit[1];
    fromMatrix[7] = fromUnit[0];
    fromMatrix[8] = 0;


    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 1, 3, 1, fromMatrix, 3, toUnit, 1, 0.0, w, 1);

    mkl_free(fromMatrix);

    double sin = cblas_dnrm2(3, w, 1);

    double *Finv = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );
    Finv[0] = fromUnit[0];
    Finv[1] = fromUnit[1];
    Finv[2] = fromUnit[2];
    Finv[3] = v[0];
    Finv[4] = v[1];
    Finv[5] = v[2];
    Finv[6] = w[0];
    Finv[7] = w[1];
    Finv[8] = w[2];

    mkl_free(w);
    mkl_free(v);

    double *G = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );
    G[0] = dot;
    G[1] = -sin;
    G[2] = 0;
    G[3] = sin;
    G[4] = dot;
    G[5] = 0;
    G[6] = 0;
    G[7] = 0;
    G[8] = 1;

    double *scrap = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, Finv, 3, G, 3, 0, scrap, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                3, 3, 3, 1, scrap, 3, Finv, 3, 0, out, 3);

    mkl_free(scrap);
    mkl_free(G);
    mkl_free(Finv);
};


