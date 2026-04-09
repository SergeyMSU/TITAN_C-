#pragma once
#include "Header.h"
using namespace std;


class CoolingFunction
{
public:
    // Данные функции охлаждения
    std::vector<double> T_values;          // температура [K]
    std::vector<double> Lambda_values;     // Lambda/n_H^2 [erg cm^3/s]
    std::vector<double> lnT;               // логарифмы T
    std::vector<double> lnLambda;          // логарифмы Lambda

    void ReadCoolingFunction(const std::string& filename);

    double InterpolateCooling(double T) const;
    void WriteCheckFile(const string& name) const;

};

