//
// Created by ronald on 26-11-19.
//

#ifndef LAPACKIK_LINK_H
#define LAPACKIK_LINK_H

class Link {

protected:
    double jointAngle;
    double radius;
    double twist;
    double offset;

public:
    Link(double jointAngle, double radius, double offset, double twist) :
        jointAngle(jointAngle), radius(radius), offset(offset), twist(twist){
    }

    double* getForwardMatrix();
};


#endif //LAPACKIK_LINK_H