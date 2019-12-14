//
// Created by ronald on 14-12-19.
//

#ifndef LAPACKIK_SCHUNK_H
#define LAPACKIK_SCHUNK_H

#include <vector>
#include "Link.h"


class Schunk {

private:
    std::vector<Link*> *links = new std::vector<Link*>(7);

public:
    Schunk(std::vector<Link*> *links): links(links){};
    void getForwardMatrix(double* out);
};


#endif //LAPACKIK_SCHUNK_H
