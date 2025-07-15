#pragma once
#include "Header.h"


struct PrintOptions {
    std::optional<double> x;
    std::optional<double> y;
    std::optional<double> z;
    std::optional<string> fluid;
};

class Phys_param
{
public:

    vector<string> param_names;
    // ��� ����� ����������, ������� ������� �� ������� � �������!
    vector<string> H_name;   // ����� �������������� ��������� ��������
    vector<string> pui_name;   // ����� �������������� ��������� �������
    vector<string> MK_param;   // �������������� ��������� ��� MK

    bool is_PUI = false;       // ������� �� ������?


    double Velosity_inf = -2.54327;   // �������� �������� ����� �� �������������
    double B_inf = 8.91006;           // �������� ������ ���������� ���� �� �������������
    double alphaB_inf = 0.76604444;     // ����������� ���������� ���� �� �������������
    double B_0 = 114.037;             // ��������� ���� �� 1 �.�.
    double p_0 = 3118.94;             // �������� �� 1 �.�.
    double par_a_2 = 0.102046;         // �������� � ������� �����������
    double par_n_H_LISM = 3.0;        // �/� ������������ �������� �� ������������� 
    double par_Kn = 39.3412;
    double mn_He_0 = 0.035194;        // ��������� ������������ ����� �� 1 ��
    double mn_He_inf = 0.15;        // ��������� ������������ ����� �� inf
    short int num_H = 4;           // ������� ������ �������� ��� ��

    // ���������
    double gamma = (5.0 / 3.0);       // ���������� ��������
    double g1 = (5.0 / 3.0 - 1.0);       // ���������� ��������

    // ����������� ���������
    double R_0 = 0.233017;             // ����������� ������ 1 �.�.
    double char_rho = 0.06;             // ����������� ������������ � ���  ��^-3
    double char_v = 10.3804;         // ����������� �������� � ��/�


    // ��������� ������� ������
    double KFL = 0.8;                   // �������� �������
    bool TVD = true;                   // ������ �� ���?

    bool culc_plasma = true;           // ������� �� ������? ����� ���������� ������ ��� ������� ��������
    bool move_HP = true;               // ������� �� HP
    bool move_BS = true;

    bool sglag_TS = true;              // ������ �� ����������� TS
    double velocity_TS = 0.01;
    double sglag_TS_k = 0.005;         // ����������� �� ������� �������
    //double sglag_TS_k_sphere = 0.001; // 0.3;   // ����������� � �������� � ��������� �����
    double sglag_TS_k_sphere_head = 0.02; // 0.08;   // ����������� � �������� �����
    double sglag_TS_k_sphere_tail = 0.01; // 0.03;   // ����������� � ��������� �����
    // �� ����� ���� ������ 1

    bool sglag_HP = true;
    double velocity_HP = 0.1;
    double sglag_HP_k_sphere = 0.002;  //0.001    // C���������� � �������� �����
    double sglag_HP_k = 0.001;          // ����������� �� � �������� �������
    double sglag_HP_angle = 1.2;    // ����������� �������� ����������� �� ����
    double sglag_HP_along = 1.0;    // ����������� �������� ����������� ����� �
    double sglag_HP_sphere = 5.0;   // ����������� �������� ����������� � �������� �������


    bool sglag_BS = false;
    double sglag_BS_k = 0.05;


    bool null_bn_on_HP = true;   // ��� ����� ����� � HP �������� ���������� ���������� ���������� ����
    bool bn_in_p_on_HP = true;   // ��� ����� ����� � HP ���������� ��������� ���� � ��������


     // ��������� ������� ��  ********************************************************
    
    bool save_AMR;        // ����� �� ��������� ����������� ������� �������������?
    bool culc_AMR;        // ����� �� ������� ������� �������������?
    bool refine_AMR;      // ����� �� �������� ����������� ������� �������������?
    unsigned int N_per_gran;  // ������� � ������� ������ �������� � ������ �����
    bool culc_cell_moments;    // ����� �� ������� ������� � �������?
    bool de_refine_AMR;        // ����� �� ��������� AMR �����, ���� ���������?
    string MK_file;            // ���� ��� ���������� �������� � �������
     

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

    void Get_Potok(const double& rho, const double& p, const double& u, const double& v,
        const double& w, const double& bx, const double& by, const double& bz,
        const double& n1, const double& n2, const double& n3, 
        std::vector<double>& Potok);  // Potok - 8-������ ������

    double energy(const double& rho, const double& p, const double& vv, const double& bb);

    void chlld(unsigned short int n_state, // �����
        const double& al, const double& be, const double& ge, // �������
        const double& w, // �������� �����
        const std::vector<double>& qqq1, const std::vector<double>& qqq2, // �������� ��������� � ���� ������ 
        std::vector<double>& qqq, // �������� �����
        bool null_bn,  // ����  true, �� ����� �������� ����� ���������� ���� ����� �����������
        unsigned short int n_disc,    // ������� ��� ����������� ��������
        const std::vector<double>& konvect_left,  // �������������� ���������� ������������� �������� �����
        const std::vector<double>& konvect_right, // �������������� ���������� ������������� �������� ������
        std::vector<double>& konvect, // �������������� ���������� ������������� �������� ������
        double& dsr, double& dsc, double& dsl,
        PrintOptions& opts, bool left_ydar = false);

    // ���� �������� ********
    void lev(const double& enI, const double& pI, const double& rI, const double& enII,//
        const double& pII, const double& rII, double& uuu, double& fee);
    void devtwo(const double& enI, const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
        const double& w, double& p);
    void newton(const double& enI, const double& pI, const double& rI, const double& enII, const double& pII, const double& rII, //
        const double& w, double& p);
    void Godunov_Solver_Alexashov(const std::vector<double>& qqq1, const std::vector<double>& qqq2,//
        const std::vector<double>& n, std::vector<double>& qqq,//
        double& dsl, double& dsp, double& dsc, double w, bool contact = false);

    // ***********

    void raspad_testing(void);
};

