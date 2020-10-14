//
// Created by ronald on 20-12-19.
//

#ifndef LAPACKIK_INVERSEOPERATIONS_H
#define LAPACKIK_INVERSEOPERATIONS_H

#include "Schunk.h"

void computeWristPosition(double *endEffectorPosition, double *endEffectorOrientation, const Schunk *schunk, double *out);

void computeWristPositionWRTShoulder(double *w, const Schunk *schunk, double *out);

double computeJointElbowAngle(const Schunk *schunk, double* sLw);

void computeSkW(double *sLw, double *out);

void computeRotationPsi(double *skw, double psi, double *out);

double computeBaseAngleShadow(double *w);

double computeShoulderAngleShadow(const Schunk *schunk, double * sLw, double *w);

void computeReferencePlaneRotation(double baseAngleSchadow, double shoulderAngleShadow, double *out);

void computeXs(double * sKw, double * referencePlaneRotation, double * out);

void computeYs(double * sKw, double * referencePlaneRotation, double * out);

void computeZs(double * sKw, double * referencePlaneRotation, double * out);

double computeBaseAngle(double psi, double * Xs, double * Ys, double * Zs);

double computeFirstShoulderAngle(double psi, double * Xs, double * Ys, double * Zs);

double computeSecondShoulderAngle(double psi, double * Xs, double * Ys, double * Zs);

void computeShoulderRotation(double * rotationMatrix, double * referencePlaneRotation, double * out);

void computeXw(double * elbowRotation, double * sKw, double * referencePlaneRotation, double * efRotation, double * out);

void computeYw(double * elbowRotation, double * sKw, double * referencePlaneRotation, double * efRotation, double * out);

void computeZw(double * elbowRotation, double * sKw, double * referencePlaneRotation, double * efRotation, double * out);

double computeFirstWristAngle(double psi, double * Xw, double * Yw, double * Zw);

double computeSecondWristAngle(double psi, double * Xw, double * Yw, double * Zw, double firstWristAngle);

double computeEFAngle(double psi, double * Xw, double * Yw, double * Zw);

void test4R7(double *Xw, double *Yw, double *Zw, double psi, double * M4RR7);

#endif //LAPACKIK_INVERSEOPERATIONS_H
