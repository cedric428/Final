//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include <iostream>
#include "BoxMuller.h"

std::vector<double> BoxMuller::getNormals(unsigned int n) {
    std::vector<double> normals = std::vector<double>();
    unsigned int i = 0;
    while (i < n) {
        double u1;
        double u2;
        double x;
        do {
            u1 = uniformGenerator.getNextUniform();
            u2 = uniformGenerator.getNextUniform();
            u1 = 2 * u1 - 1;
            u2 = 2 * u2 - 1;
            x = u1 * u1 + u2 * u2;
            i += 2;
        } while (x > 1);
        double y = sqrt(-2 * log(x) / x);
        double z1 = u1 * y;
        double z2 = u2 * y;
        if (i <= n) {
            normals.push_back(z1);
            normals.push_back(z2);
        }
    }
    return normals;
}
