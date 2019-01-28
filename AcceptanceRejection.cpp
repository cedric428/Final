//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include "AcceptanceRejection.h"

std::vector<double> AcceptanceRejection::getNormals(unsigned int n) {
    std::vector<double> normals = std::vector<double>();
    unsigned int i = 0;
    while (i < n) {
        double u1;
        double u2;
        double u3;
        double x;
        do {
            u1 = uniformGenerator.getNextUniform();
            u2 = uniformGenerator.getNextUniform();
            u3 = uniformGenerator.getNextUniform();
            i += 3;
            x = -log(u1);
        } while (u2 > exp(-(x - 1) * (x - 1) / 2));
        if (u3 <= 0.5) {
            x = -x;
        }
        if (i < n) {
            normals.push_back(x);
        }
    }
    return normals;
}
