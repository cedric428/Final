//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_INVERSETRANSFORMMETHOD_H
#define MONTECARLO_INVERSETRANSFORMMETHOD_H


#include "LinearCongruentialGenerator.h"

class InverseTransformMethod {
private:
    LinearCongruentialGenerator uniformGenerator = LinearCongruentialGenerator();
    const double a0 = 2.50662823884;
    const double a1 = -18.61500062529;
    const double a2 = 41.39119773534;
    const double a3 = -25.44106049637;

    const double b0 = -8.47351093090;
    const double b1 = 23.08336743743;
    const double b2 = -21.06224101826;
    const double b3 = 3.13082909833;

    const double c0 = 0.3374754822726147;
    const double c1 = 0.9761690190917186;
    const double c2 = 0.1607979714918209;
    const double c3 = 0.0276438810333863;
    const double c4 = 0.0038405729373609;
    const double c5 = 0.0003951896511919;
    const double c6 = 0.0000321767881768;
    const double c7 = 0.0000002888167364;
    const double c8 = 0.0000003960315187;

public:
    double getNextNormal();
};


#endif //MONTECARLO_INVERSETRANSFORMMETHOD_H
