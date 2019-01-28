//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_BOXMULLER_H
#define MONTECARLO_BOXMULLER_H


#include "LinearCongruentialGenerator.h"
#include "NormalGenerator.h"

class BoxMuller : public NormalGenerator {
private:
    LinearCongruentialGenerator uniformGenerator = LinearCongruentialGenerator();
public:
    std::vector<double> getNormals(unsigned int n) override;
};


#endif //MONTECARLO_BOXMULLER_H
