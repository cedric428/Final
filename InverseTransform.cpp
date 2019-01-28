//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include "InverseTransform.h"

double InverseTransform::getNextNormal() {
    auto u = uniformGenerator.getNextUniform();
    double y = u - 0.5;
    double x;
    if (std::abs(y) < 0.42) {
        double r = y * y;
        x = y * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1);
    } else {
        double r = u;
        if (y > 0) {
            r = 1 - u;
        }
        r = log(-log(r));
        x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
        if (y < 0) {
            x = -x;
        }
    }
    return x;
}

std::vector<double> InverseTransform::getNormals(unsigned int n) {
    std::vector<double> normals = std::vector<double>();
    unsigned int i = 0;
    while (i++ < n) {
        normals.push_back(getNextNormal());
    }
    return normals;
}
