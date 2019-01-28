//
// Created by Harish Reddy on 9/16/2017.
//

#include "LinearCongruentialGenerator.h"

double LinearCongruentialGenerator::getNextUniform() {
    x_prev = (a * x_prev + c) % k;
    return double(x_prev) / k;
}
