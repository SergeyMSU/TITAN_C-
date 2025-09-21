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
    vector<string> plasma_name; // ����� ���������� ��������
    vector<string> plasma_pui_name; // ����� ���������� ��������
    vector<string> H_param_names;


    vector<string> pui_name;   // ����� �������������� ��������� �������
    vector<string> MK_param;   // �������������� ��������� ��� MK

    vector<string> p_pui_name; 
    vector<string> H_name;   // ����� �������������� ��������� ��������

    std::unordered_set<std::string> r2_snos_names;
    std::unordered_set<std::string> r2g_snos_names;

    unordered_map<string, double> perevod_razmer;


    bool is_PUI;       // ������� �� ������?

    uint8_t num_pui;           // ������� ������ ������� � �������

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>  pui_in_zone;   
    // ������� ���������� � ����� ���� ����� ����� ������� ������������


    // ��������� ������� ������ ��� ������� ���� ����������� (�� ���)
    // Inner_Hard, Outer_Hard, Outer_Soft
    // 0 - ���� ������� �� �����, ��� ��� ��� ��� ���� ����������
    // 1 - ��������   cells[0] | cells[0]
    // 2 - ������  cells[0] | gran
    // 3 - ����������� ������ cells[0] | Center 

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  pui_condition;
    // ������� ��������� ������� ��� �������

    Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>  plasma_condition;
    // ������� ��������� ������� ��� ������

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_condition;
    // ������� ��������� ������� ��� ��������

    // ������� �������������� ������ �� ��������
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_1;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_1;

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_2;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_2;

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_3;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_3;

    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  hydrogen_arise_4;
    Eigen::Matrix< int8_t, Eigen::Dynamic, Eigen::Dynamic>  proton_arise_4;

    // �������� ������ proton_arise � ���, ��� � ��� ��� "����������� ��������"
    // �� ���� ���� � ������� 1 ��� pui_2, �� � ������� ��� ������ �������
    // ��� ������� ������, ��� ��� pui_2 - ��� ��������� �������
    // �� ���� � ������� �� ����� pui, ������� �� � ��������� �������, �� ���� ������������
    // �������� ���������� ����������.

    std::unordered_map<std::string, bool> Culc_hidrogen;   
    // ����� ���� �������� ���������, � ����� ���
    // ��� ���������� ���� �� ������� ������ ������� �����-�� ��������


    bool is_div_V_in_cell;       // ������� �� ����������� �������� � ������� 
    // �� ������ ����? �����, ��������, ��� ����� �������
    // ����� � ��� ������� ������������, ������� ��� ��������� ������



    std::function<void (const short int&,
        unordered_map<string, double>&,
        unordered_map<string, double>&)> Plasma_components;



    double Velosity_inf = -2.54327;   // �������� �������� ����� �� �������������
    double B_inf = 8.91006;           // �������� ������ ���������� ���� �� �������������
    double alphaB_inf = 0.6981317;     // ���� � �������� - ����������� ���������� ���� �� �������������
    double B_0 = 114.037;             // ��������� ���� �� 1 �.�.
    double p_0 = 3118.94;             // �������� �� 1 �.�.
    double par_a_2 = 0.102046;         // �������� � ������� �����������
    double par_n_H_LISM = 3.0;        // �/� ������������ �������� �� ������������� 
    double par_Kn = 40.0906;
    double mrho_He_0 = 0.035194;        // ��������� ������������ ����� �� 1 ��
    double mrho_He_inf = 0.15;        // ��������� ������������ ����� �� inf
    double mep = 0.000544617;      // ��������� ����� ��������� � ����� �������
    uint8_t num_H = 4;           // ������� ������ �������� ��� ��

    // ���������
    double gamma = (5.0 / 3.0);       // ���������� ��������
    double g1 = (5.0 / 3.0 - 1.0);       // ���������� �������� - 1

    // ����������� ���������
    double R_0 = 0.237455;             // ����������� ������ 1 �.�.
    double char_rho = 0.06;             // ����������� ������������ � ���  ��^-3
    double char_v = 10.3804;         // ����������� �������� � ��/�


    // ��������� ������� ������
    double KFL;                   // �������� �������
    bool TVD;                   // ������ �� ��� ��� ������?
    bool TVD_atom = true;                   // ������ �� ��� ��� ������?

    bool Snos_on_HP = false;       // ������ �� �������� ���� �� HP � ���� ������ (����� ������� ������ ��������)

    bool culc_plasma;     // �� �������      // ������� �� ������? ����� ���������� ������ ��� ������� ��������
    bool culc_atoms;             // ��������� �� ����� ��� ��������� �� ������������
    bool move_setka;
    
    
    bool move_TS;               // ������� �� HP
    bool move_HP;               // ������� �� HP
    bool move_BS;

    bool sglag_TS;              // ������ �� ����������� TS
    double velocity_TS;
    double sglag_TS_k;         // ����������� �� ������� �������
    double sglag_TS_k_sphere;         // ����������� �� ������� �������
    //double sglag_TS_k_sphere = 0.001; // 0.3;   // ����������� � �������� � ��������� �����
    double sglag_TS_k_sphere_head; // 0.08;   // ����������� � �������� �����
    double sglag_TS_k_sphere_tail; // 0.03;   // ����������� � ��������� �����
    // �� ����� ���� ������ 1

    bool sglag_HP;
    double velocity_HP;
    double sglag_HP_k_sphere;  //0.005 0.002    // C���������� � �������� �����
    double sglag_HP_k_angle;  //0.005 0.002    // C���������� � ��������� ����� � ��������� �����������
    double sglag_HP_k; // 0.001         // ����������� �� � �������� �������
    double sglag_HP_angle;    // 1.2 ����������� �������� ����������� �� ����
    double sglag_HP_along;    // ����������� �������� ����������� ����� �
    double sglag_HP_sphere;   // ����������� �������� ����������� � �������� ������� - �� �������


    bool sglag_BS;
    double sglag_BS_k;
    double velocity_BS;


    bool null_bn_on_HP;   // ��� ����� ����� � HP �������� ���������� ���������� ���������� ����
    bool bn_in_p_on_HP;   // ��� ����� ����� � HP ���������� ��������� ���� � ��������
    bool contact_hard;    // ���� �� �� �������� ������������ �������� ����� � �������� ��������
    bool TS_hard;    // ���� �� �� �������� ������������ �������� ����� � �������� ��������

     // ��������� ������� ��  ********************************************************
    
    bool save_AMR;        // ����� �� ��������� ����������� ������� �������������?
    bool culc_AMR;        // ����� �� ������� ������� �������������?
    bool refine_AMR;      // ����� �� �������� ����������� ������� �������������?
    unsigned int N_per_gran;  // ������� � ������� ������ �������� � ������ �����
    bool culc_cell_moments;    // ����� �� ������� ������� � �������?
    bool de_refine_AMR;        // ����� �� ��������� AMR �����, ���� ���������?
    string MK_file;            // ���� ��� ���������� �������� � �������

    bool MK_source_S = false;         // ����� �� ������� ��������� S+ � S- � ��
    int pui_nW;                  // ����� ��������� �������� S+ � S- � �������
    double pui_wR;               // ������������ ������������� �������� � �������� S+ � S-
    short int pui_h0_n = 300;    // ��� h0 � ��������� pui
    double pui_h0_wc = 200;      // ��� h0 � ��������� pui
    short int pui_F_n = 100;      // �� ������� ������ �� ��������� ������������� ��� ��������� PUI
     

    //inline double sigma(double x);                  // ������� �����������
    //inline double sigma2(double x, double y);

    inline double sigma(double x)
    {
        return kv(1.0 - this->par_a_2 * log(x));
    }

    inline double sigma2(double x, double y)
    {
        return kv(1.0 - this->par_a_2 * log((x) * (y)));
    }


    Eigen::Matrix3d Matr;              // ������� �������� 1
    Eigen::Matrix3d Matr2;             // ������� �������� 2


    // ������ ��� ���������� ����� ���������� �� ��������� �����
    std::vector<double> heliolat_deg;
    std::vector<double> n_p_cm3;
    std::vector<double> V_kms;
    std::vector<double> T_K;


    Phys_param();
    void set_parameters(void);

    double Get_rho_0(const double& the);
    double Get_v_0(const double& the);
    double Get_T_0(const double& the);

    double Get_razmer(string par);  // ���������� ����������� ��� �������� 
    // ��������� � ��������� ��������

    // ���������� ��������� ������
    void Plasma_components_1(const short int& zone,
        unordered_map<string, double>& param_in_cell,
        unordered_map<string, double>& param);

    void Plasma_components_2(const short int& zone,
        unordered_map<string, double>& param_in_cell,
        unordered_map<string, double>& param);


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
        PrintOptions& opts, bool left_ydar = false, bool contact = false);

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

private:
    using VarRef = std::variant<
        double*,
        int*,
        uint8_t*,
        unsigned int*,
        bool*,
        std::string*
    >;

    std::unordered_map<std::string, VarRef> varMap;

    void initVarMap();
    void parseAndAssign(VarRef varRef, const std::string& valueStr);
};

