//
// Created by ronald on 20-12-19.
//

#include <mkl.h>
#include "Schunk.h"
#include "InverseOperations.h"


void computeWristPosition(double *endEffectorPosition, double *endEffectorOrientation, const Schunk *schunk, double *out){

    double *zUnit = (double *)mkl_malloc( 3*1*sizeof( double ), 64 );
    zUnit[0] = 0;
    zUnit[1] = 0;
    zUnit[2] = 1;

    double *scrap = (double *)mkl_malloc( 3*1*sizeof( double ), 64 );

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 1, 3, schunk->getWlt(), endEffectorOrientation, 3, zUnit, 1, 0.0, scrap, 1);

    for(int i=0; i<3; i++){
        out[i] = endEffectorPosition[i] - scrap[i];
    }

    mkl_free(scrap);
    mkl_free(zUnit);
}

void computeWristPositionWRTShoulder(double *w, const Schunk *schunk, double *out){

    out[0] = w[0];
    out[1] = w[1];
    out[2] = w[2] - schunk->getBls();
}

double computeJointElbowAngle(const Schunk *schunk, double* sLw){

    double norm = cblas_dnrm2(3, sLw, 1);

    double thetaPrime = acos((pow(schunk->getSle(),2) + pow(schunk->getElw(),2) - pow(norm,2))/(2*schunk->getSle()*schunk->getElw()));
    return M_PI - thetaPrime;
}

void computeSkW(double *sLw, double *out){

    double norm = cblas_dnrm2(3, sLw, 1);
    double *suw = (double *)mkl_malloc( 3*1*sizeof( double ), 64 );

    for(int i=0;i<3;i++){
        suw[i] = sLw[i] / norm;
    }

    out[0] = 0;
    out[1] = -suw[2];
    out[2] = suw[1];
    out[3] = suw[2];
    out[4] = 0;
    out[5] = -suw[0];
    out[6] = -suw[1];
    out[7] = suw[0];
    out[8] = 0;
}

void computeRotationPsi(double *skw, double psi, double *out){

    double *eye = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );

    eye[0] = 1;
    eye[1] = 0;
    eye[2] = 0;
    eye[3] = 0;
    eye[4] = 1;
    eye[5] = 0;
    eye[6] = 0;
    eye[7] = 0;
    eye[8] = 1;

    double *scrap = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 1, 3, 1, skw, 3, skw, 3, 0.0, scrap, 3);

    for(int i=0;i<3*3;i++){
        out[i] = eye[i] + (1-cos(psi)) * scrap[i] + sin(psi) * skw[i];
    }

    mkl_free(scrap);
    mkl_free(eye);
}

double computeBaseAngleShadow(double *w){
    return atan2(w[1], w[0]);
}

double computeShoulderAngleShadow(const Schunk *schunk, double * sLw, double *w){

    double euclidSlw = cblas_dnrm2(3, sLw, 1);

    double alpha = asin((w[2]- schunk->getBls())/euclidSlw);
    double beta = acos((pow(schunk->getSle(), 2) + pow(euclidSlw, 2) - pow(schunk->getElw(), 2))/(2*schunk->getSle()*euclidSlw));

    return M_PI_2 - alpha - beta;
}

void computeReferencePlaneRotation(double baseAngleSchadow, double shoulderAngleShadow, double *out){

    out[0] = cos(baseAngleSchadow) * cos(shoulderAngleShadow);
    out[1] = -cos(baseAngleSchadow) * sin(shoulderAngleShadow);
    out[2] = -sin(baseAngleSchadow);
    out[3] = cos(shoulderAngleShadow) * sin(baseAngleSchadow);
    out[4] = -sin(baseAngleSchadow) * sin(shoulderAngleShadow);
    out[5] = cos(baseAngleSchadow);
    out[6] = -sin(shoulderAngleShadow);
    out[7] = -cos(shoulderAngleShadow);
    out[8] = 0;
}

void computeXs(double * sKw, double * referencePlaneRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, sKw, 3, referencePlaneRotation, 3, 0, out, 3);
}

void computeYs(double * sKw, double * referencePlaneRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, sKw, 3, sKw, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, -1, out, 3, referencePlaneRotation, 3, 0, out, 3);
}

void computeZs(double * sKw, double * referencePlaneRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, sKw, 3, sKw, 3, 0, out, 3);

    out[0]++;
    out[4]++;
    out[8]++;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, referencePlaneRotation, 3, 0, out, 3);
}

double computeBaseAngle(double psi, double * Xs, double * Ys, double * Zs){

    return atan((-sin(psi)*Xs[4] -cos(psi)*Ys[4] - Zs[4])/(-sin(psi)*Xs[1] -cos(psi)*Ys[1] - Zs[1]));
}

double computeFirstShoulderAngle(double psi, double * Xs, double * Ys, double * Zs){

    return acos(-sin(psi)*Xs[7] -cos(psi)*Ys[7] -Zs[7]);
}

double computeSecondShoulderAngle(double psi, double * Xs, double * Ys, double * Zs){

    return atan((sin(psi)*Xs[8] + cos(psi)*Ys[8] + Zs[8])/(-sin(psi)*Xs[6] -cos(psi)*Ys[6] - Zs[6]));
}

void computeShoulderRotation(double * rotationMatrix, double * referencePlaneRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, rotationMatrix, 3, referencePlaneRotation, 3, 0, out, 3);
}

void computeXw(double * elbowRotation, double * sKw, double * referencePlaneRotation, double * efRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            3, 3, 3, 1, sKw, 3, referencePlaneRotation, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, elbowRotation, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, efRotation, 3, 0, out, 3);

}

void computeYw(double * elbowRotation, double * sKw, double * referencePlaneRotation, double * efRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, sKw, 3, sKw, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, referencePlaneRotation, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, elbowRotation, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                3, 3, 3, -1, out, 3, efRotation, 3, 0, out, 3);
}

void computeZw(double * elbowRotation, double * sKw, double * referencePlaneRotation, double * efRotation, double * out){

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, sKw, 3, sKw, 3, 0, out, 3);

    out[0] += 1;
    out[4] += 1;
    out[8] += 1;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, referencePlaneRotation, 3, 0, out, 3);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, elbowRotation, 3, 0, out, 3);


    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                3, 3, 3, 1, out, 3, efRotation, 3, 0, out, 3);
}

double computeFirstWristAngle(double psi, double * Xw, double * Yw, double * Zw, double secondWristAngle){

    if(secondWristAngle < 0){
        return atan2(-(sin(psi)*Xw[5] + cos(psi)*Yw[5] +Zw[5]), -(sin(psi)*Xw[2] + cos(psi)*Yw[2] + Zw[2]))
    }else{
        return atan2((sin(psi)*Xw[5] + cos(psi)*Yw[5] +Zw[5]), (sin(psi)*Xw[2] + cos(psi)*Yw[2] + Zw[2]));
    }
}

double computeSecondWristAngle(double psi, double * Xw, double * Yw, double * Zw){

    return -acos(sin(psi)*Xw[8] + cos(psi)*Yw[8] + Zw);
}

double computeEFAngle(double psi, double * Xw, double * Yw, double * Zw, double secondWristAngle){

    if(secondWristAngle < 0){
        return atan2(-(sin(psi)*Xw[7] + cos(psi)*Yw[7] + Zw[7]),(sin(psi)*Xw[6] + cos(psi)*Yw[6] + Zw[6]));
    }else{
        return atan2((sin(psi)*Xw[7] + cos(psi)*Yw[7] + Zw[7]),-(sin(psi)*Xw[6] + cos(psi)*Yw[6] + Zw[6]));
    }
}

void test4R7(double *Xw, double *Yw, double *Zw, double psi, double *M4RR7) {

    double* res = (double *)mkl_malloc( 3*3*sizeof( double ), 64 );

    for(int i = 0;i<9;i++){

        res[i] = sin(psi) * Xw[i] + cos(psi) * Yw[i] + Zw[i] - M4RR7[i];
    }

    mkl_free(res);

}
