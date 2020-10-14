//
// Created by ronald on 26-11-19.
//

#ifndef LAPACKIK_ROTATIONJOINT_H
#define LAPACKIK_ROTATIONJOINT_H


#include "Link.h"

class RotationJoint : public Link{

public:
    RotationJoint(double jointAngle, double radius, double offset, double twist) :
            Link(jointAngle, radius, offset, twist){}
    void setJointAngle(double jointAngle);
};


#endif //LAPACKIK_ROTATIONJOINT_H
