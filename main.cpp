/* C source code is found in dgemm_example.c */

#include <stdio.h>
#include <iostream>
#include <vector>
#include "mkl.h"
#include "model/Link.h"
#include "model/Schunk.h"
#include <math.h>
#include "matplotlibcpp.h"
#include "model/operations.h"
#include "model/InverseOperations.h"


int main(){

    double psi = 0;

    double bls = 0.3;
    double sle = 0.328;
    double elw = 0.323;
    double wlt = 0.0824;

    double initAngle1 = 0;
    double initAngle2 = M_PI_4;
    double initAngle3 = 0;
    double initAngle4 = M_PI_2;
    double initAngle5 = 0;
    double initAngle6 = -M_PI_4;
    double initAngle7 = 0;

    std::vector<double> *initAngles = new std::vector<double>();

    initAngles->push_back(initAngle1);
    initAngles->push_back(initAngle2);
    initAngles->push_back(initAngle3);
    initAngles->push_back(initAngle4);
    initAngles->push_back(initAngle5);
    initAngles->push_back(initAngle6);
    initAngles->push_back(initAngle7);

    Schunk *schunk = new Schunk(bls, sle, elw, wlt);
    schunk->setJointAngles(*initAngles);

    double *forwardMatrix = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    schunk->getForwardMatrix(forwardMatrix);

    double *efPosition = (double *)mkl_malloc( 3*1*sizeof( double ), 64 );
    double *efOrientation = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );

    efPosition[0] = forwardMatrix[3];
    efPosition[1] = forwardMatrix[7];
    efPosition[2] = forwardMatrix[11];

    efOrientation[0] = forwardMatrix[0];
    efOrientation[1] = forwardMatrix[1];
    efOrientation[2] = forwardMatrix[2];
    efOrientation[3] = forwardMatrix[4];
    efOrientation[4] = forwardMatrix[5];
    efOrientation[5] = forwardMatrix[6];
    efOrientation[6] = forwardMatrix[8];
    efOrientation[7] = forwardMatrix[9];
    efOrientation[8] = forwardMatrix[10];

    mkl_free(forwardMatrix);

    double *wristPosition = (double *)mkl_malloc( 3*1*sizeof(double), 64);
    computeWristPosition(efPosition, efOrientation, schunk, wristPosition);

    double *sLw = (double*)mkl_malloc(3*1* sizeof(double), 64);
    computeWristPositionWRTShoulder(wristPosition, schunk, sLw);

    double elbowAngle = computeJointElbowAngle(schunk, sLw);

    double *skw = (double*)mkl_malloc(3*3* sizeof(double), 64);
    computeSkW(sLw, skw);

    double *rotationPsi = (double*)mkl_malloc(3*3*(sizeof(double)), 64);
    computeRotationPsi(skw, psi, rotationPsi);

    double baseAngleShadow = computeBaseAngleShadow(wristPosition);
    double shoulderAngleShadow = computeShoulderAngleShadow(schunk, sLw, wristPosition);

    double *referencePlaneRotation = (double*)mkl_malloc(3*3* sizeof(double), 64);
    computeReferencePlaneRotation(baseAngleShadow, shoulderAngleShadow, referencePlaneRotation);

    double *Xs = (double*)mkl_malloc(3*3* sizeof(double), 64);
    double *Ys = (double*)mkl_malloc(3*3* sizeof(double), 64);
    double *Zs = (double*)mkl_malloc(3*3* sizeof(double), 64);

    computeXs(skw, referencePlaneRotation, Xs);
    computeYs(skw, referencePlaneRotation, Ys);
    computeZs(skw, referencePlaneRotation, Zs);

    double baseAngle = computeBaseAngle(psi, Xs, Ys, Zs);
    double firstShoulderAngle = computeFirstShoulderAngle(psi, Xs, Ys, Zs);
    double secondShoulderAngle = computeSecondShoulderAngle(psi, Xs, Ys, Zs);

    double * shoulderRotation = (double *)mkl_malloc(3*3*sizeof(double), 64);
    computeShoulderRotation(rotationPsi, referencePlaneRotation, shoulderRotation);

    double * Xw = (double*)mkl_malloc(3*3*sizeof(double), 64);
    double * Yw = (double*)mkl_malloc(3*3*sizeof(double), 64);
    double * Zw = (double*)mkl_malloc(3*3*sizeof(double), 64);

    double * threeRFour = (double*)mkl_malloc(3*3* sizeof(double), 64);
    schunk->getLinks()->at(3)->getRotation(threeRFour);

    double * M4R7 = (double*)mkl_malloc(4*43*3*sizeof(double), 64);
    schunk->get4R7(M4R7);

    computeXw(threeRFour, skw, referencePlaneRotation, efOrientation, Xw);
    computeYw(threeRFour, skw, referencePlaneRotation, efOrientation, Yw);
    computeZw(threeRFour, skw, referencePlaneRotation, efOrientation, Zw);

    test4R7(Xw, Yw, Zw, psi, M4R7);

    double firstWristAngle = computeFirstWristAngle(psi, Xw, Yw, Zw);
    double secondWristAngle = computeSecondWristAngle(psi, Xw, Yw, Zw, firstWristAngle);
    double efAngle = computeEFAngle(psi, Xw, Yw, Zw);

    std::vector<double>* resultingJointAngles = new std::vector<double>(7);
    resultingJointAngles->push_back(baseAngle);
    resultingJointAngles->push_back(firstShoulderAngle);
    resultingJointAngles->push_back(secondShoulderAngle);
    resultingJointAngles->push_back(elbowAngle);
    resultingJointAngles->push_back(firstWristAngle);
    resultingJointAngles->push_back(secondWristAngle);
    resultingJointAngles->push_back(efAngle);

    Schunk *schunk2 = new Schunk(bls, sle, elw, wlt);
    schunk2->setJointAngles(*resultingJointAngles);

    double *forwardMatrix2 = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );
    schunk2->getForwardMatrix(forwardMatrix2);

    return 0;
}
