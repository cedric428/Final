

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "BlackScholes.h"
#include "MonteCarloNonPathDependant.h"
#include "MonteCarloPathDependant.h"
#include "AcceptanceRejection.h"
#include "BoxMuller.h"

void runMonteCarloForPlainVanillaEuropeanOptions();

void runMonteCarloForPathDependantOptions();

void runMonteCarloToCompareRNGs();

int main() {
    // Problem 1
    runMonteCarloForPlainVanillaEuropeanOptions();

    // Problem 2
    runMonteCarloForPathDependantOptions();

    // Problem 3
    runMonteCarloToCompareRNGs();

    std::cout << "The output files are : p1.csv, p2.csv, p3.csv" << std::endl;
    return 0;
}

void runMonteCarloForPlainVanillaEuropeanOptions() {
    std::cout << "Running problem 1." << std::endl;
    BlackScholes blackScholes = BlackScholes(41, 42, 0.03, 0.01, 0.2, 0.75);
    auto bsCallValue = blackScholes.getCallValue();
    auto bsCallDelta = blackScholes.getCallDelta();
    auto bsCallVega = blackScholes.getCallVega();
    auto bsPutValue = blackScholes.getPutValue();
    auto bsPutDelta = blackScholes.getPutDelta();
    auto bsPutVega = blackScholes.getPutVega();

    MonteCarloNonPathDependant monteCarlo = MonteCarloNonPathDependant(new InverseTransform(), 41, 42, 0.03, 0.01, 0.2,
                                                                       0.75);
    monteCarlo.runSimulation();
    unsigned int k;
    std::ofstream outFile;
    outFile.open("p1.csv");
    for (k = 0; k < 10; k++) {
        unsigned long long int N = 10000 * pow(2, k);
        auto mcCallValue = monteCarlo.getCallValue(N);
        auto mcCallDelta = monteCarlo.getCallDelta(N);
        auto mcCallVega = monteCarlo.getCallVega(N);
        outFile << std::setprecision(14) << N << ","
                << mcCallValue << "," << sqrt(N) * fabs(bsCallValue - mcCallValue) << ","
                << mcCallDelta << "," << sqrt(N) * fabs(bsCallDelta - mcCallDelta) << ","
                << mcCallVega << "," << sqrt(N) * fabs(bsCallVega - mcCallVega)
                << std::endl;
    }
    outFile << std::endl;
    for (k = 0; k < 10; k++) {
        unsigned long long int N = 10000 * pow(2, k);
        auto mcPutValue = monteCarlo.getPutValue(N);
        auto mcPutDelta = monteCarlo.getPutDelta(N);
        auto mcPutVega = monteCarlo.getPutVega(N);
        outFile << std::setprecision(14) << N << ","
                << mcPutValue << "," << sqrt(N) * fabs(bsPutValue - mcPutValue) << ","
                << mcPutDelta << "," << sqrt(N) * fabs(bsPutDelta - mcPutDelta) << ","
                << mcPutVega << "," << sqrt(N) * fabs(bsPutVega - mcPutVega)
                << std::endl;
    }
    outFile.close();
    std::cout << "Finished problem 1." << std::endl << std::endl;
}

void runMonteCarloForPathDependantOptions() {
    std::cout << "Running problem 2." << std::endl;
    BlackScholes blackScholes = BlackScholes(39, 39, 0.02, 0.01, 0.25, 0.75);
    auto bsValue = blackScholes.getDAOValue(35);

    MonteCarloPathDependant monteCarloFixedTimesteps = MonteCarloPathDependant(new InverseTransform(), 39, 39, 0.02,
                                                                               0.01, 0.25, 0.75, 35);
    monteCarloFixedTimesteps.runSimulation(200);
    unsigned int k;
    std::ofstream outFile;
    outFile.open("p2.csv");
    for (k = 0; k < 10; k++) {
        unsigned long long int N = 10000 * pow(2, k);
        unsigned long long int mk = ceil(pow(N, 1.0 / 3) * pow(0.75, 2.0 / 3));
        unsigned long long int nk = floor(N / mk);
        MonteCarloPathDependant monteCarloVariableTimesteps = MonteCarloPathDependant(new InverseTransform(), 39, 39,
                                                                                      0.02, 0.01, 0.25, 0.75, 35);
        monteCarloVariableTimesteps.runSimulation(mk);
        auto mcValueVariableTimestep = monteCarloVariableTimesteps.getDAOPrice(N);
        auto mcValueFixedTimestep = monteCarloFixedTimesteps.getDAOPrice(N);
        outFile << std::setprecision(14) << N << "," << 200 << "," << N / 200 << ","
                << mcValueFixedTimestep << "," << fabs(bsValue - mcValueFixedTimestep) << ","
                << N << "," << mk << "," << nk << ","
                << mcValueVariableTimestep << "," << fabs(bsValue - mcValueVariableTimestep)
                << std::endl;
    }
    outFile.close();
    std::cout << "Finished problem 2." << std::endl << std::endl;
}

void runMonteCarloToCompareRNGs() {
    std::cout << "Running problem 3." << std::endl;
    BlackScholes blackScholes = BlackScholes(50, 55, 0.04, 0.0, 0.3, 0.5);
    auto bsValue = blackScholes.getPutValue();

    MonteCarloNonPathDependant monteCarloInverseTransform = MonteCarloNonPathDependant(new InverseTransform(), 50, 55,
                                                                                       0.04, 0.0, 0.3, 0.5);
    monteCarloInverseTransform.runSimulation();
    MonteCarloNonPathDependant monteCarloAcceptanceRejection = MonteCarloNonPathDependant(new AcceptanceRejection(), 50,
                                                                                          55, 0.04, 0.0, 0.3, 0.5);
    monteCarloAcceptanceRejection.runSimulation();
    MonteCarloNonPathDependant monteCarloBoxMuller = MonteCarloNonPathDependant(new BoxMuller(), 50, 55, 0.04, 0.0, 0.3,
                                                                                0.5);
    monteCarloBoxMuller.runSimulation();

    unsigned int k;
    std::ofstream outFile;
    outFile.open("p3.csv");
    for (k = 0; k < 10; k++) {
        unsigned long long int N = 10000 * pow(2, k);
        InverseTransform temp1 = InverseTransform();
        AcceptanceRejection temp2 = AcceptanceRejection();
        BoxMuller temp3 = BoxMuller();

        // Calculate how many normals are generated using N values from uniform distribution for each method.
        auto mcValueInverseTransform = monteCarloInverseTransform.getPutValue(temp1.getNormals(N).size());
        auto nar = temp2.getNormals(N).size();
        auto mcValueAcceptanceRejection = monteCarloAcceptanceRejection.getPutValue(nar);
        auto nbm = temp3.getNormals(N).size();
        auto mcValueBoxMuller = monteCarloBoxMuller.getPutValue(nbm);

        outFile << std::setprecision(14) << N << ","
                << mcValueInverseTransform << "," << fabs(bsValue - mcValueInverseTransform) << ","
                << nar << "," << mcValueAcceptanceRejection << "," << fabs(bsValue - mcValueAcceptanceRejection) << ","
                << nbm << "," << mcValueBoxMuller << "," << fabs(bsValue - mcValueBoxMuller)
                << std::endl;

    }
    outFile.close();
    std::cout << "Finished problem 3." << std::endl << std::endl;
}
