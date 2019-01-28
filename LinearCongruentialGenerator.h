//
// Created by Harish Reddy on 9/16/2017.
//

#ifndef MONTECARLO_LINEARCONGRUENTIALGENERATOR_H
#define MONTECARLO_LINEARCONGRUENTIALGENERATOR_H


class LinearCongruentialGenerator {
public:
    double getNextUniform();

private:
    const unsigned int a = 39373;
    const unsigned int c = 0;
    const unsigned int k = (1 << 31) - 1;
    unsigned long long int x_prev = 1;
};


#endif //MONTECARLO_LINEARCONGRUENTIALGENERATOR_H
