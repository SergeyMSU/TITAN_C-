#include <filesystem>
#include <string>
#include <utility>
#include <cstdint>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>
#include <iterator>
#include <cstdlib>
#include <optional>
#include <omp.h>
#include <chrono>
#include <random>

using namespace std;

#define kv(x) ((x) * (x))

const double const_pi = 3.14159265358979323846;

// Задаём сечение перезарядки
double sig(const double& u)
{
    const double a_2 = 0.1307345665;
    return  (kv(1.0 - a_2 * log(u)));
}

void Culc_nu_maxwel(void)
{
    double U, cp;

    double UL = 0.01;
    double UR = 200.0;
    unsigned int UN = 1000;

    double cpL = 0.01;
    double cpR = 100.0;
    unsigned int cpN = 100;

    double R, R0, al;
    unsigned int N1 = 100;
    unsigned int N2 = 100;
    double S, dR0, dal, SS, SSS;

    std::ofstream outfile("Int_1.bin", std::ios::binary);
    std::ofstream outfile2("Int_1.txt");

    outfile2 << "TITLE = HP  VARIABLES = U, cp, I, I_malama, err" << endl;

    for (int i = 0; i < UN; ++i)
    {       
        for (int j = 0; j < cpN; ++j) 
        {  
            // Тело вложенного цикла
            U = UL + i * (UR - UL) / (UN - 1);
            cp = cpL + j * (cpR - cpL) / (cpN - 1);

            S = 0.0;
            SS = 0.0;
            SSS = 0.0;
            dal = const_pi / N2;
            dR0 = 5.0 * cp / N1;

#pragma omp parallel for reduction(+:S, SS, SSS) private(R0, al, R)
            for (int ii = 0; ii < N1; ++ii)
            {
                R0 = (ii + 0.5) * 5.0 * cp / N1;

                for (int jj = 0; jj < N2; ++jj)
                {
                    al = (jj + 0.5) * const_pi / N2;
                    R = sqrt(kv(R0 * sin(al)) + kv(U - R0 * cos(al)));
                    S += R * sig(R) * kv(R0) * sin(al) * exp(-kv(R0 / cp));
                    SS += R * kv(R0) * sin(al) * exp(-kv(R0 / cp));
                    SSS += kv(R0) * sin(al) * exp(-kv(R0 / cp));
                }
            }

            S *= 2.0 * const_pi * dal * dR0;
            SS *= 2.0 * const_pi * dal * dR0;
            SSS *= 2.0 * const_pi * dal * dR0;
            S /= pow(sqrt(const_pi) * cp, 3);
            SS /= pow(sqrt(const_pi) * cp, 3);
            SSS /= pow(sqrt(const_pi) * cp, 3);
            SS *= sig(SS/SSS);
            outfile.write(reinterpret_cast<const char*>(&S), sizeof(S));
            outfile2 << U << " " << cp << " " << S << " " << SS << " " << SS * 100/S - 100.0 << endl;
        }
    }

    outfile.close();
    outfile2.close();
}


int main()
{
    Culc_nu_maxwel();
    std::cout << "Hello World!\n";
}
