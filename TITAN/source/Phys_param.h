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
    // Все имена параметров, наличие которых мы ожидаем в ячейках!
    vector<string> H_name;   // имена дополнительных жидкостей водорода
    vector<string> pui_name;   // имена дополнительных жидкостей пикапов
    vector<string> MK_param;   // дополнительные параметры для MK

    bool is_PUI = false;       // Считаем ли пикапы?


    double Velosity_inf = -2.54327;   // Значение скорости смеси на бесконечности
    double B_inf = 8.91006;           // Значение МОДУЛЯ магнитного поля на бесконечности
    double alphaB_inf = 0.76604444;     // Направление магнитного поля на бесконечности
    double B_0 = 114.037;             // Магнитное поле на 1 а.е.
    double p_0 = 3118.94;             // давление на 1 а.е.
    double par_a_2 = 0.102046;         // Параметр в сечении перезарядки
    double par_n_H_LISM = 3.0;        // б/р Концентрация водорода на бесконечности 
    double par_Kn = 39.3412;
    double mn_He_0 = 0.035194;        // Множитель концентрации гелия на 1 АЕ
    double mn_He_inf = 0.15;        // Множитель концентрации гелия на inf
    short int num_H = 4;           // Сколько сортов водорода для МК

    // параметры
    double gamma = (5.0 / 3.0);       // показатель адиабаты
    double g1 = (5.0 / 3.0 - 1.0);       // показатель адиабаты

    // Характерные параметры
    double R_0 = 0.233017;             // Характерный размер 1 а.е.
    double char_rho = 0.06;             // Характерная концентрация в СГС  см^-3
    double char_v = 10.3804;         // Характерная скорость в км/с


    // Настройки расчёта Плазмы
    double KFL = 0.8;                   // критерий Куранта
    bool TVD = true;                   // Делаем ли ТВД?

    bool culc_plasma = true;           // Считаем ли плазму? Можно заморозить плазму для расчёта водорода
    bool move_HP = true;               // Двигаем ли HP
    bool move_BS = true;

    bool sglag_TS = true;              // Делаем ли сглаживание TS
    double velocity_TS = 0.01;
    double sglag_TS_k = 0.005;         // Сглаживание на высоких широтах
    //double sglag_TS_k_sphere = 0.001; // 0.3;   // Сглаживание в головной и хвостовой части
    double sglag_TS_k_sphere_head = 0.02; // 0.08;   // Сглаживание в головной части
    double sglag_TS_k_sphere_tail = 0.01; // 0.03;   // Сглаживание в хвостовой части
    // Не может быть больше 1

    bool sglag_HP = true;
    double velocity_HP = 0.1;
    double sglag_HP_k_sphere = 0.002;  //0.001    // Cглиживание в головной части
    double sglag_HP_k = 0.001;          // Сглаживание не в головной области
    double sglag_HP_angle = 1.2;    // коэффициент усилинея сглаживания по углу
    double sglag_HP_along = 1.0;    // коэффициент усилинея сглаживания вдоль х
    double sglag_HP_sphere = 5.0;   // коэффициент усиления сглаживания в головной области


    bool sglag_BS = false;
    double sglag_BS_k = 0.05;


    bool null_bn_on_HP = true;   // Для ячеек рядом с HP обнуляем нормальную компоненту магнитного поля
    bool bn_in_p_on_HP = true;   // Для ячеек рядом с HP записываем магнитное поле в давление


     // Настройки расчёта МК  ********************************************************
    
    bool save_AMR;        // Нужно ли сохранять посчитанные функции распределения?
    bool culc_AMR;        // Нужно ли считать функции распределения?
    bool refine_AMR;      // Нужно ли мельчить посчитанные функции распределения?
    unsigned int N_per_gran;  // Сколько в среднем частиц вылетает с каждой грани
    bool culc_cell_moments;    // Нужно ли считать моменты в ячейках?
    bool de_refine_AMR;        // Нужно ли огрублять AMR сетку, если требуется?
    string MK_file;            // Фаил для сохранения моментов в ячейках
     

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

    void Get_Potok(const double& rho, const double& p, const double& u, const double& v,
        const double& w, const double& bx, const double& by, const double& bz,
        const double& n1, const double& n2, const double& n3, 
        std::vector<double>& Potok);  // Potok - 8-мерный вектор

    double energy(const double& rho, const double& p, const double& vv, const double& bb);

    void chlld(unsigned short int n_state, // метод
        const double& al, const double& be, const double& ge, // нормаль
        const double& w, // скорость грани
        const std::vector<double>& qqq1, const std::vector<double>& qqq2, // основные параметры с двух сторон 
        std::vector<double>& qqq, // Выходной поток
        bool null_bn,  // если  true, то нужно обнулить поток магнитного поля через поверхность
        unsigned short int n_disc,    // формула для определения скорости
        const std::vector<double>& konvect_left,  // Дополнительные переменные конвективного переноса слева
        const std::vector<double>& konvect_right, // Дополнительные переменные конвективного переноса справа
        std::vector<double>& konvect, // Дополнительные переменные конвективного переноса ПОТОКИ
        double& dsr, double& dsc, double& dsl,
        PrintOptions& opts, bool left_ydar = false);

    // Блок Годунова ********
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

