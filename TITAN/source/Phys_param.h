#pragma once
#include "Header.h"

class Phys_param
{
public:
    double Velosity_inf = -2.54327;   // �������� �������� ����� �� �������������
    double B_inf = 8.91006;           // �������� ������ ���������� ���� �� �������������
    double alphaB_inf = 0.698132;     // ����������� ���������� ���� �� �������������
    double B_0 = 114.037;             // ��������� ���� �� 1 �.�.
    double p_0 = 3118.94;             // �������� �� 1 �.�.

    // ����������� ���������
    double R_0 = 0.233017;             // ����������� ������ 1 �.�.
    double char_rho = 0.06;             // ����������� ������������ � ���  ��^-3
    double char_v = 10.3804;         // ����������� �������� � ��/�


    Eigen::Matrix3d Matr;              // ������� �������� 1
    Eigen::Matrix3d Matr2;             // ������� �������� 2




    // ������ ��� ���������� ����� ���������� �� ��������� �����
    std::vector<double> heliolat_deg;
    std::vector<double> n_p_cm3;
    std::vector<double> V_kms;
    std::vector<double> T_K;


    Phys_param();

    double Get_rho_0(const double& the);
    double Get_v_0(const double& the);
};

