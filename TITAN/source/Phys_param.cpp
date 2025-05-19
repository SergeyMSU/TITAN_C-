#include "Phys_param.h"

Phys_param::Phys_param()
{
	this->Matr << -0.9958639688067077,  0.01776569097515556,  0.08910295088675518,
		           0.07561695085992419, 0.7057402284561812,   0.7044237408557894,
		          -0.05036896241933166, 0.7082479157489926,  -0.7041646522383864;

	this->Matr2 << -0.9958639688067080, 0.0756169508599243, -0.0503689624193315,
	            	0.0177656909751554, 0.7057402284561816,  0.7082479157489927,
		            0.0891029508867553, 0.7044237408557898, -0.7041646522383865;

    // Открываем файл для чтения
    std::ifstream file("nVT1au_new_2000-2022av.dat");
    if (!file.is_open()) {
        std::cerr << "Error 6564397608" << std::endl;
        exit(-1);
    }

    std::string line;
    bool is_first_line = true;

    // Читаем файл построчно
    while (std::getline(file, line)) {
        // Пропускаем первую строку (заголовок)
        if (is_first_line) {
            is_first_line = false;
            continue;
        }

        // Разбиваем строку на отдельные значения
        std::istringstream iss(line);
        double hlat, np, v, t;

        // Если строка содержит данные, записываем их в векторы
        if (iss >> hlat >> np >> v >> t) {
            this->heliolat_deg.push_back(hlat);
            this->n_p_cm3.push_back(np);
            this->V_kms.push_back(v);
            this->T_K.push_back(t);
        }
    }

    file.close();
}

double Phys_param::Get_rho_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->n_p_cm3;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1)) / this->char_rho;
}

double Phys_param::Get_v_0(const double& the)
{
    const auto& angles = this->heliolat_deg;
    const auto& np = this->V_kms;

    // Проверка на выход за пределы диапазона
    if (the < angles.front() || the > angles.back()) {
        cout << "Ygol " + std::to_string(the) + " out of diapazon!" << endl;
        cout << "ERROR  0989767567" << endl;
    }

    // Поиск ближайших точек
    size_t i = 0;
    while (i < angles.size() - 1 && angles[i + 1] < the) {
        i++;
    }

    // Если угол точно совпадает с одним из значений в файле
    if (std::abs(the - angles[i]) < 1e-6) {
        return np[i];
    }

    // Линейная интерполяция: n_p = np1 + (theta - theta1) * (np2 - np1) / (theta2 - theta1)
    double theta1 = angles[i];
    double theta2 = angles[i + 1];
    double np1 = np[i];
    double np2 = np[i + 1];

    return (np1 + (the - theta1) * (np2 - np1) / (theta2 - theta1)) / this->char_v;

    return 0.0;
}

void Phys_param::Get_Potok(const double& rho, const double& p, const double& u, const double& v, 
    const double& w, const double& bx, const double& by, const double& bz, 
    const double& n1, const double& n2, const double& n3, std::vector<double>& Potok)
{
    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d m;

    n << n1, n2, n3;

    if (Potok.size() < 8 || std::abs(n.norm() - 1.0) > 1e-6)
    {
        cout << "Error 8765086751   " << n.norm() << "  " << Potok.size() << endl;
        exit(-1);
    }

    get_bazis(n, t, m);


    double vn = scalarProductFast(u, v, w, n1, n2, n3);
    double vt = scalarProductFast(u, v, w, t(0), t(1), t(2));
    double vm = scalarProductFast(u, v, w, m(0), m(1), m(2));

    double Bn = scalarProductFast(bx, by, bz, n1, n2, n3);
    double Bt = scalarProductFast(bx, by, bz, t(0), t(1), t(2));
    double Bm = scalarProductFast(bx, by, bz, m(0), m(1), m(2));

    double bb = kvv(bx, by, bz);
    double vv = kvv(u, v, w);

    double e = this->energy(rho, p, vv, bb);
    double pT = p + bb / 2.0;

    Potok[0] = rho * vn; 
    Potok[7] = (e + pT) * vn - Bn * scalarProductFast(bx, by, bz, u, v, w);
    double Pvn = rho * kv(vn) + pT - kv(Bn);
    double Pvt = rho * vn * vt - Bn * Bt;
    double Pvm = rho * vn * vm - Bn * Bm;
    double Pbt = Bt * vn - Bn * vt;
    double Pbm = Bm * vn - Bn * vm;
    Potok[1] = Pvn * n(0) + Pvt * t(0) + Pvm * m(0);
    Potok[2] = Pvn * n(1) + Pvt * t(1) + Pvm * m(1);
    Potok[3] = Pvn * n(2) + Pvt * t(2) + Pvm * m(2);
    Potok[4] = Pbt * t(0) + Pbm * m(0);
    Potok[5] = Pbt * t(1) + Pbm * m(1);
    Potok[6] = Pbt * t(2) + Pbm * m(2);
}

double Phys_param::energy(const double& rho, const double& p, const double& vv, const double& bb)
{
    return p / (this->gamma - 1) + rho * vv / 2.0 + bb / 2.0;
}


void Phys_param::chlld(unsigned short int n_state, // метод
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
    PrintOptions& opts)  // Дополнительные опциональные параметры
    // n_state = 1 HLL, // 2 HLLC,  3 HLLD
{

    double FL[8];
    double FR[8];
    double UL[8];
    double UR[8];
    double UZ[8];

    double wv = w;

    double r1 = qqq1[0];
    double u1 = qqq1[1];
    double v1 = qqq1[2];
    double w1 = qqq1[3];
    double p1 = qqq1[4];
    double bx1 = qqq1[5] / spi4;
    double by1 = qqq1[6] / spi4;
    double bz1 = qqq1[7] / spi4;

    double r2 = qqq2[0];
    double u2 = qqq2[1];
    double v2 = qqq2[2];
    double w2 = qqq2[3];
    double p2 = qqq2[4];
    double bx2 = qqq2[5] / spi4;
    double by2 = qqq2[6] / spi4;
    double bz2 = qqq2[7] / spi4;

    double ro = (r2 + r1) / 2.0;
    double au = (u2 + u1) / 2.0;
    double av = (v2 + v1) / 2.0;
    double aw = (w2 + w1) / 2.0;
    double ap = (p2 + p1) / 2.0;
    double abx = (bx2 + bx1) / 2.0;
    double aby = (by2 + by1) / 2.0;
    double abz = (bz2 + bz1) / 2.0;

    double bk = abx * al + aby * be + abz * ge;
    double b2 = kv(abx) + kv(aby) + kv(abz);

    Eigen::Vector3d n;
    Eigen::Vector3d t;
    Eigen::Vector3d m;

    n << al, be, ge;

    get_bazis(n, t, m);

    Eigen::Vector3d vL;
    Eigen::Vector3d vR;
    Eigen::Vector3d bL;
    Eigen::Vector3d bR;


    vL(0) = scalarProductFast(n(0), n(1), n(2), u1, v1, w1);
    vL(1) = scalarProductFast(t(0), t(1), t(2), u1, v1, w1);
    vL(2) = scalarProductFast(m(0), m(1), m(2), u1, v1, w1);
    
    vR(0) = scalarProductFast(n(0), n(1), n(2), u2, v2, w2);
    vR(1) = scalarProductFast(t(0), t(1), t(2), u2, v2, w2);
    vR(2) = scalarProductFast(m(0), m(1), m(2), u2, v2, w2);

    bL(0) = scalarProductFast(n(0), n(1), n(2), bx1, by1, bz1);
    bL(1) = scalarProductFast(t(0), t(1), t(2), bx1, by1, bz1);
    bL(2) = scalarProductFast(m(0), m(1), m(2), bx1, by1, bz1);

    bR(0) = scalarProductFast(n(0), n(1), n(2), bx2, by2, bz2);
    bR(1) = scalarProductFast(t(0), t(1), t(2), bx2, by2, bz2);
    bR(2) = scalarProductFast(m(0), m(1), m(2), bx2, by2, bz2);


    double aaL = bL(0) / sqrt(r1);
    double b2L = bL.dot(bL);
    double b21 = b2L / r1;
    double cL = sqrt(this->gamma * p1 / r1);
    double qp = sqrt(b21 + cL * (cL + 2.0 * aaL));
    double qm = sqrt(b21 + cL * (cL - 2.0 * aaL));
    double cfL = (qp + qm) / 2.0;
    double ptL = p1 + b2L / 2.0;

    double aaR = bR(0) / sqrt(r2);
    double b2R = bR.dot(bR);
    double b22 = b2R / r2;
    double cR = sqrt(this->gamma * p2 / r2);
    qp = sqrt(b22 + cR * (cR + 2.0 * aaR));
    qm = sqrt(b22 + cR * (cR - 2.0 * aaR));
    double cfR = (qp + qm) / 2.0;
    double ptR = p2 + b2R / 2.0;

    double aC = (aaL + aaR) / 2.0;
    double b2o = (b22 + b21) / 2.0;
    double cC = sqrt(this->gamma * ap / ro);
    qp = sqrt(b2o + cC * (cC + 2.0 * aC));
    qm = sqrt(b2o + cC * (cC - 2.0 * aC));
    double cfC = (qp + qm) / 2.0;
    double vC1 = (vL(0) + vR(0)) / 2.0;

    double SL = 0.0, SR = 0.0;

    if (n_disc == 0)
    {
        SL = min((vL(0) - cfL), (vR(0) - cfR));
        SR = max((vL(0) + cfL), (vR(0) + cfR));
    }
    else if (n_disc == 1)
    {
        SL = min((vL(0) - cfL), (vC1 - cfC));
        SR = max((vR(0) + cfR), (vC1 + cfC));
    }
    else if (n_disc == 2)
    {
        double SL_1 = min((vL(0) - cfL), (vC1 - cfC));
        double SR_1 = max((vR(0) + cfR), (vC1 + cfC));
        double SL_2 = min((vL(0) - cfL), (vR(0) - cfR));
        double SR_2 = max((vL(0) + cfL), (vR(0) + cfR));
        double oo = 0.6;
        double oo1 = 1.0 - oo;
        SL = oo * SL_1 + oo1 * SL_2;
        SR = oo * SR_1 + oo1 * SR_2;
    }
    else
    {
        cout << "Error 8654354789" << endl;
    }

    dsr = SR;
    dsl = SL;

    double suR = SR - vR(0);
    double suL = SL - vL(0);
    double SM = (suR * r2 * vR(0) - ptR + ptL - suL * r1 * vL(0))
        / (suR * r2 - suL * r1);

    double upt1 = (kv(u1) + kv(v1) + kv(w1)) / 2.0;
    double sbv1 = u1 * bx1 + v1 * by1 + w1 * bz1;

    double upt2 = (kv(u2) + kv(v2) + kv(w2)) / 2.0;
    double sbv2 = u2 * bx2 + v2 * by2 + w2 * bz2;

    double e1 = p1 / g1 + r1 * upt1 + b2L / 2.0;
    double e2 = p2 / g1 + r2 * upt2 + b2R / 2.0;

    FL[0] = r1 * vL[0];
    FL[1] = r1 * vL[0] * vL[0] + ptL - kv(bL[0]);
    FL[2] = r1 * vL[0] * vL[1] - bL[0] * bL[1];
    FL[3] = r1 * vL[0] * vL[2] - bL[0] * bL[2];
    FL[4] = (e1 + ptL) * vL[0] - bL[0] * sbv1;
    FL[5] = 0.0;
    FL[6] = vL[0] * bL[1] - vL[1] * bL[0];
    FL[7] = vL[0] * bL[2] - vL[2] * bL[0];

    FR[0] = r2 * vR[0];
    FR[1] = r2 * vR[0] * vR[0] + ptR - kv(bR[0]);
    FR[2] = r2 * vR[0] * vR[1] - bR[0] * bR[1];
    FR[3] = r2 * vR[0] * vR[2] - bR[0] * bR[2];
    FR[4] = (e2 + ptR) * vR[0] - bR[0] * sbv2;
    FR[5] = 0.0;
    FR[6] = vR[0] * bR[1] - vR[1] * bR[0];
    FR[7] = vR[0] * bR[2] - vR[2] * bR[0];

    UL[0] = r1;
    UL[4] = e1;
    UR[0] = r2;
    UR[4] = e2;

    for (short unsigned int i = 0; i < 3; i++)
    {
        UL[i + 1] = r1 * vL[i];
        UL[i + 5] = bL[i];
        UR[i + 1] = r2 * vR[i];
        UR[i + 5] = bR[i];
    }

    for (short unsigned int ik = 0; ik < 8; ik++)
    {
        UZ[ik] = (SR * UR[ik] - SL * UL[ik] + FL[ik] - FR[ik]) / (SR - SL);
    }

    if (null_bn == true) UZ[5] = 0.0;

    if (n_state == 1)
    {
        double dq[8];
        double FW[8];

        for (short unsigned int ik = 0; ik < 8; ik++)
        {
            dq[ik] = UR[ik] - UL[ik];
        }

        double TL = SL;
        double TR = SR;
        if (SL > wv)
        {
            TL = 0.0;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UL[ik];
            }
        }
        else if (SL <= wv && wv <= SR)
        {
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UZ[ik];
            }
        }
        else if (SR < wv)
        {
            TR = 0.0;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UR[ik];
            }
        }
        else
        {
            cout << "error  8767854390  " << endl;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                cout << "ik: " << qqq1[ik] << "   " << qqq2[ik] << ";" << endl;
            }
            exit(-1);
        }

        double a = TR * TL;
        double b = TR - TL;

        Eigen::Vector3d qv;
        Eigen::Vector3d qb;

        qqq[0] = (TR * FL[0] - TL * FR[0] + a * dq[0]) / b - FW[0];
        qqq[8] = (TR * FL[8] - TL * FR[8] + a * dq[8]) / b - FW[8];
        qqq[4] = (TR * FL[4] - TL * FR[4] + a * dq[4]) / b - FW[4];
        for (short unsigned int ik = 1; ik < 4; ik++)
        {
            qv(ik - 1) = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }
        for (short unsigned int ik = 5; ik < 8; ik++)
        {
            qb(ik - 5) = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }

        for (short unsigned int ik = 0; ik < 3; ik++)
        {
            qqq[ik + 1] = n(ik) * qv(0) + t(ik) * qv(1) + m(ik) * qv(2);
            qqq[ik + 5] = n(ik) * qb(0) + t(ik) * qb(1) + m(ik) * qb(2);
        }

        return;
    }

    return;
}