//
// Created by Harish Reddy on 9/17/2017.
//

#include <cmath>
#include "AcceptanceRejectionMethod.h"

double AcceptanceRejectionMethod::getNextNormal() {
    double u1;
    double u2;
    double u3;
    double x;
    do {
        u1 = uniformGenerator.getNextUniform();
        u2 = uniformGenerator.getNextUniform();
        u3 = uniformGenerator.getNextUniform();
        x = -log(u1);
    } while (u2 > exp(-(x - 1) * (x - 1) / 2));
    if (u3 <= 0.5) {
        x = -x;
    }
    return x;
}
