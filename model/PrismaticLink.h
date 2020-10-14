//
// Created by ronald on 26-11-19.
//

#ifndef LAPACKIK_PRISMATICLINK_H
#define LAPACKIK_PRISMATICLINK_H


#include "Link.h"

class PrismaticLink: public Link {

public:
    PrismaticLink(double jointAngle, double radius, double offset, double twist) :
            Link(jointAngle, radius, offset, twist){}
    void setTranslationRadius(double radius);
};


#endif //LAPACKIK_PRISMATICLINK_H
