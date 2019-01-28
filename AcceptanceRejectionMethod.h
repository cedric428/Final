//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_ACCEPTANCEREJECTIONMETHOD_H
#define MONTECARLO_ACCEPTANCEREJECTIONMETHOD_H


#include "LinearCongruentialGenerator.h"

class AcceptanceRejectionMethod {
private:
    LinearCongruentialGenerator uniformGenerator = LinearCongruentialGenerator();

public:
    double getNextNormal();
};


#endif //MONTECARLO_ACCEPTANCEREJECTIONMETHOD_H
