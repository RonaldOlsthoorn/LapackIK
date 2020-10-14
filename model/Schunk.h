//
// Created by ronald on 14-12-19.
//

#ifndef LAPACKIK_SCHUNK_H
#define LAPACKIK_SCHUNK_H

#include <vector>
#include <math.h>
#include "Link.h"


class Schunk {

private:
    double bls;
    double sle;
    double elw;
    double wlt;

    std::vector<Link*> *links;

public:
    Schunk(double bls, double sle, double elw, double wlt): bls(bls), sle(sle), elw(elw), wlt(wlt){

        Link* l1 = new Link(0, 0, bls, -M_PI/2);
        Link* l2 = new Link(0, 0, 0, M_PI/2);
        Link* l3 = new Link(0, 0, sle, -M_PI/2);
        Link* l4 = new Link(0, 0, 0, M_PI/2);
        Link* l5 = new Link(0, 0, elw, -M_PI/2);
        Link* l6 = new Link(0, 0, 0, M_PI/2);
        Link* l7 = new Link(0, 0, wlt, 0);

        links = new std::vector<Link*>();

        links->push_back(l1);
        links->push_back(l2);
        links->push_back(l3);
        links->push_back(l4);
        links->push_back(l5);
        links->push_back(l6);
        links->push_back(l7);
    };

    void getForwardMatrix(double* out) const;
    void setJointAngles(std::vector<double> angles);
    std::vector<double> getEndEffectorPosition();

    double getBls() const {return bls;}
    double getSle() const {return sle;}
    double getElw() const {return elw;}
    double getWlt() const {return wlt;}

    std::vector<Link*> * getLinks() const {return links;}

    void get4R7(double* out) const;

};


#endif //LAPACKIK_SCHUNK_H
