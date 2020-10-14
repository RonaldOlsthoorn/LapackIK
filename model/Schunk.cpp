//
// Created by ronald on 14-12-19.
//

#include <mkl_service.h>
#include <mkl.h>
#include <cstdio>
#include "Schunk.h"

void Schunk::getForwardMatrix(double *out) const {

    links->front()->getForwardMatrix(out);

    for(std::vector<Link*>::iterator it = ++(links->begin()); it != links->end(); it++){
        double* forwardMatrix = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );
        (*it)->getForwardMatrix(forwardMatrix);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    4, 4, 4, 1, out, 4, forwardMatrix, 4, 0, out, 4);

        mkl_free(forwardMatrix);
    }
}

void Schunk::setJointAngles(std::vector<double> angles) {

    int i = 0;

    for(auto &link : *links){
        link->setJointAngle(angles.at(i));
        i++;
    }

}

std::vector<double> Schunk::getEndEffectorPosition() {

    double* forwardMatrix = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    getForwardMatrix(forwardMatrix);

    double* origin = (double *)mkl_malloc( 4*1*sizeof( double ), 64 );

    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;
    origin[3] = 1;

    double* out = (double *)mkl_malloc( 4*1*sizeof( double ), 64 );

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                4, 1, 4, 1.0, forwardMatrix, 4, origin, 1, 0.0, out, 1);


    mkl_free(forwardMatrix);
    mkl_free(origin);

    std::vector<double> res(3);

    res.push_back(out[0]);
    res.push_back(out[1]);
    res.push_back(out[2]);

    mkl_free(out);

    return res;
}

void Schunk::get4R7(double *out) const {

    double* m5 = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );
    double* m6 = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );
    double* m7 = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );

    links->at(4)->getRotation(m5);
    links->at(5)->getRotation(m6);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, m5, 3, m6, 3, 0, out, 3);

    links->at(6)->getRotation(m7);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, m7, 3, 0, out, 3);

    mkl_free(m5);
    mkl_free(m6);
    mkl_free(m7);
}
