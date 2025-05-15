#pragma once
#include "Header.h"

class Phys_param
{
public:
    double Velosity_inf = -2.54327;   // Значение скорости смеси на бесконечности
    double B_inf = 8.91006;           // Значение МОДУЛЯ магнитного поля на бесконечности
    double alphaB_inf = 0.698132;     // Направление магнитного поля на бесконечности
    double B_0 = 114.037;             // Магнитное поле на 1 а.е.
    double p_0 = 3118.94;             // давление на 1 а.е.

    // Характерные параметры
    double R_0 = 0.233017;             // Характерный размер 1 а.е.
    double char_rho = 0.06;             // Характерная концентрация в СГС  см^-3
    double char_v = 10.3804;         // Характерная скорость в км/с


    Eigen::Matrix3d Matr;              // Матрица перехода 1
    Eigen::Matrix3d Matr2;             // Матрица перехода 2




    // Данные для считывания файла параметров на начальной сфере
    std::vector<double> heliolat_deg;
    std::vector<double> n_p_cm3;
    std::vector<double> V_kms;
    std::vector<double> T_K;


    Phys_param();

    double Get_rho_0(const double& the);
    double Get_v_0(const double& the);
};

