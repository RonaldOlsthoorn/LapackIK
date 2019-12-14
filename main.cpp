/* C source code is found in dgemm_example.c */

#include <stdio.h>
#include <iostream>
#include <vector>
#include "mkl.h"
#include "model/Link.h"
#include "model/Schunk.h"
#include <math.h>


int min(int x, int y){
    if(x<y){
        return x;
    }
    return y;
}

int main(){

    double bls = 0.3;
    double sle = 0.328;
    double elw = 0.323;
    double wlt = 0.0824;

    std::vector<Link*> *links = new std::vector<Link*>();

    Link* l1 = new Link(0, 0, bls, -M_PI/2);
    Link* l2 = new Link(0, 0, 0, M_PI/2);
    Link* l3 = new Link(0, 0, sle, -M_PI/2);
    Link* l4 = new Link(0, 0, 0, M_PI/2);
    Link* l5 = new Link(0, 0, elw, -M_PI/2);
    Link* l6 = new Link(0, 0, 0, M_PI/2);
    Link* l7 = new Link(0, 0, wlt, 0);


    links->push_back(l1);
    links->push_back(l2);
    links->push_back(l3);
    links->push_back(l4);
    links->push_back(l5);
    links->push_back(l6);
    links->push_back(l7);

    Schunk *schunk = new Schunk(links);

    double *runningMatrix = (double *)mkl_malloc( 4*4*sizeof( double ), 64 );

    schunk->getForwardMatrix(runningMatrix);

    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            printf ("%12.5G", runningMatrix[4*i+j]);
        }
    printf("\n");
    }


    return 0;
}
