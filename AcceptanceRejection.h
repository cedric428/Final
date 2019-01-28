//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_ACCEPTANCEREJECTIONMETHOD_H
#define MONTECARLO_ACCEPTANCEREJECTIONMETHOD_H


#include <vector>
#include "LinearCongruentialGenerator.h"
#include "NormalGenerator.h"

class AcceptanceRejection : public NormalGenerator {
private:
    LinearCongruentialGenerator uniformGenerator = LinearCongruentialGenerator();

public:
    std::vector<double> getNormals(unsigned int n) override;
};


#endif //MONTECARLO_ACCEPTANCEREJECTIONMETHOD_H
