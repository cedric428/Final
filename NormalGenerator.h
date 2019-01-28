//
// Created by Harish Reddy on 9/17/2017.
//

#ifndef MONTECARLO_NORMALGENERATOR_H
#define MONTECARLO_NORMALGENERATOR_H


#include <vector>

class NormalGenerator {
public:
    virtual std::vector<double> getNormals(unsigned int n)= 0;
};


#endif //MONTECARLO_NORMALGENERATOR_H
