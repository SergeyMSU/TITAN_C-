#include "Phys_param.h"

Phys_param::Phys_param()
{
    this->param_names.push_back("rho"); 
    this->param_names.push_back("p");
    this->param_names.push_back("Vx");
    this->param_names.push_back("Vy");
    this->param_names.push_back("Vz");
    this->param_names.push_back("Bx");
    this->param_names.push_back("By");
    this->param_names.push_back("Bz");
    this->param_names.push_back("Q");
    this->param_names.push_back("n_He");

    this->param_names.push_back("rho_H1");
    this->param_names.push_back("Vx_H1");
    this->param_names.push_back("Vy_H1");
    this->param_names.push_back("Vz_H1");
    this->param_names.push_back("p_H1");

    this->param_names.push_back("rho_H2");
    this->param_names.push_back("Vx_H2");
    this->param_names.push_back("Vy_H2");
    this->param_names.push_back("Vz_H2");
    this->param_names.push_back("p_H2");

    this->param_names.push_back("rho_H3");
    this->param_names.push_back("Vx_H3");
    this->param_names.push_back("Vy_H3");
    this->param_names.push_back("Vz_H3");
    this->param_names.push_back("p_H3");

    this->param_names.push_back("rho_H4");
    this->param_names.push_back("Vx_H4");
    this->param_names.push_back("Vy_H4");
    this->param_names.push_back("Vz_H4");
    this->param_names.push_back("p_H4");

    // Задаём имена дополнительных жидкостей
    this->H_name.push_back("_H1");
    this->H_name.push_back("_H2");
    this->H_name.push_back("_H3");
    this->H_name.push_back("_H4");



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

    double Bn = scalarProductFast(bx, by, bz, n1, n2, n3) / spi4;
    double Bt = scalarProductFast(bx, by, bz, t(0), t(1), t(2)) / spi4;
    double Bm = scalarProductFast(bx, by, bz, m(0), m(1), m(2)) / spi4;

    double bb = kvv(bx, by, bz) / spi4 / spi4;
    double vv = kvv(u, v, w);

    double e = this->energy(rho, p, vv, bb);
    double pT = p + bb / 2.0;

    Potok[0] = rho * vn; 
    Potok[7] = (e + pT) * vn - Bn * scalarProductFast(bx, by, bz, u, v, w) / spi4;
    double Pvn = rho * kv(vn) + pT - kv(Bn);
    double Pvt = rho * vn * vt - Bn * Bt;
    double Pvm = rho * vn * vm - Bn * Bm;
    double Pbt = Bt * vn - Bn * vt;
    double Pbm = Bm * vn - Bn * vm;
    Pbt *= spi4;
    Pbm *= spi4;
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
    // Конвективные переменные добавил только для HLL
{

    double FL[8];
    double FR[8];
    double UL[8];
    double UR[8];
    double UZ[8];

    vector<double> FLk(konvect_left.size());
    vector<double> FRk(konvect_left.size());
    vector<double> ULk(konvect_left.size());
    vector<double> URk(konvect_left.size());
    vector<double> UZk(konvect_left.size());

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

    //t << 0.0, 0.707107, 0.707107;
    //m << 0.0, 0.707107, -0.707107;

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
    else if (n_disc == 3)
    {
        SL = min(vL[0], vR[0]) - max(cfL, cfR);
        SR = max(vL[0], vR[0]) + max(cfL, cfR);
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
    dsc = SM;

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
    for (short int i = 0; i < konvect_left.size(); i++)
    {
        FLk[i] = konvect_left[i] * vL[0];
    }

    FR[0] = r2 * vR[0];
    FR[1] = r2 * vR[0] * vR[0] + ptR - kv(bR[0]);
    FR[2] = r2 * vR[0] * vR[1] - bR[0] * bR[1];
    FR[3] = r2 * vR[0] * vR[2] - bR[0] * bR[2];
    FR[4] = (e2 + ptR) * vR[0] - bR[0] * sbv2;
    FR[5] = 0.0;
    FR[6] = vR[0] * bR[1] - vR[1] * bR[0];
    FR[7] = vR[0] * bR[2] - vR[2] * bR[0];   // ***
    for (short int i = 0; i < konvect_right.size(); i++)
    {
        FRk[i] = konvect_right[i] * vR[0];
    }

    UL[0] = r1;
    UL[4] = e1;
    UR[0] = r2;
    UR[4] = e2;

    for (short int i = 0; i < konvect_right.size(); i++)
    {
        ULk[i] = konvect_left[i];
        URk[i] = konvect_right[i];
    }

    //whach(SL);
    //whach(SR);
    //whach(SM);


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

    for (short int ik = 0; ik < konvect_right.size(); ik++)
    {
        UZk[ik] = (SR * URk[ik] - SL * ULk[ik] + FLk[ik] - FRk[ik]) / (SR - SL);
    }

    if (null_bn == true) UZ[5] = 0.0;

    if (n_state == 1)   // HLL
    {
        double dq[8];
        double FW[8];

        vector<double> dqk(konvect_left.size());
        vector<double> FWk(konvect_left.size());

        for (short unsigned int ik = 0; ik < 8; ik++)
        {
            dq[ik] = UR[ik] - UL[ik];
        }
        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            dqk[ik] = URk[ik] - ULk[ik];
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

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                FWk[ik] = wv * ULk[ik];
            }
        }
        else if (SL <= wv && wv <= SR)
        {
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UZ[ik];
            }

            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                FWk[ik] = wv * UZk[ik];
            }
        }
        else if (SR < wv)
        {
            TR = 0.0;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                FW[ik] = wv * UR[ik];
            }
            for (short int ik = 0; ik < konvect_right.size(); ik++)
            {
                FWk[ik] = wv * URk[ik];
            }
        }
        else
        {
            cout << "error  8767854390  " << endl;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                cout << "ik: " << qqq1[ik] << "   " << qqq2[ik] << ";" << endl;
            }
            if (opts.x) whach(*opts.x);
            if (opts.y) whach(*opts.y);
            if (opts.z) whach(*opts.z);
            exit(-1);
        }

        double a = TR * TL;
        double b = TR - TL;

        Eigen::Vector3d qv;
        Eigen::Vector3d qb;

        qqq[0] = (TR * FL[0] - TL * FR[0] + a * dq[0]) / b - FW[0];
        qqq[4] = (TR * FL[4] - TL * FR[4] + a * dq[4]) / b - FW[4];
        for (short int ik = 0; ik < konvect_right.size(); ik++)
        {
            konvect[ik] = (TR * FLk[ik] - TL * FRk[ik] + a * dqk[ik]) / b - FWk[ik];
        }


        for (short unsigned int ik = 1; ik < 4; ik++)
        {
            qv(ik - 1) = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }
        for (short unsigned int ik = 5; ik < 8; ik++)
        {
            qb(ik - 5) = (TR * FL[ik] - TL * FR[ik] + a * dq[ik]) / b - FW[ik];
        }

        double SN = max(fabs(SL), fabs(SR));
        double wbn = 0.0;
        if (wv >= SR)
        {
            wbn = wv * bR[0];
        }
        else if (wv <= SL)
        {
            wbn = wv * bL[0];
        }
        else
        {
            wbn = wv * (bL[0] + bR[0]) / 2.0;
        }

        qb[0] = -SN * (bR[0] - bL[0]) - wbn;
        if (null_bn == true) qb[0] = 0.0;

        //whach(qb(0));
        //whach(qb(1));
        //whach(qb(2));

        for (short unsigned int ik = 0; ik < 3; ik++)
        {
            qqq[ik + 1] = n(ik) * qv(0) + t(ik) * qv(1) + m(ik) * qv(2);
            qqq[ik + 5] = n(ik) * qb(0) + t(ik) * qb(1) + m(ik) * qb(2);
            qqq[ik + 5] = spi4 * qqq[ik + 5];
        }

        return;
    }

    if (n_state == 2) // HLLC
    {
        double suRm = suR / (SR - SM);
        double suLm = suL / (SL - SM);
        double rzR = r2 * suRm;
        double rzL = r1 * suLm;

        Eigen::Vector3d vzR, vzL, bzR, bzL;

        vzR[0] = SM;
        vzL[0] = SM;

        double ptz = (suR * r2 * ptL - suL * r1 * ptR +
            r1 * r2 * suR * suL * (vR[0] - vL[0]))
            / (suR * r2 - suL * r1);

        bzR[0] = UZ[5];   // bR[0]
        bzL[0] = UZ[5];   // bL[0]

        vzR[1] = UZ[2] / UZ[0];
        vzR[2] = UZ[3] / UZ[0];
        vzL[1] = vzR[1];
        vzL[2] = vzR[2];

        vzR[1] = vR[1] + UZ[5] * (bR[1] - UZ[6]) / suR / r2;
        vzR[2] = vR[2] + UZ[5] * (bR[2] - UZ[7]) / suR / r2;
        vzL[1] = vL[1] + UZ[5] * (bL[1] - UZ[6]) / suL / r1;
        vzL[2] = vL[2] + UZ[5] * (bL[2] - UZ[7]) / suL / r1;

        bzR[1] = UZ[6];
        bzR[2] = UZ[7];
        bzL[1] = bzR[1];
        bzL[2] = bzR[2];

        double sbvz = (UZ[5] * UZ[1] + UZ[6] * UZ[2] + UZ[7] * UZ[3]) / UZ[0];

        double ezR = e2 * suRm + (ptz * SM - ptR * vR[0] + UZ[5] * (sbv2 - sbvz)) / (SR - SM);
        double ezL = e1 * suLm + (ptz * SM - ptL * vL[0] + UZ[5] * (sbv1 - sbvz)) / (SL - SM);

        if (fabs(UZ[5]) <= epsb)
        {
            vzR[1] = vR[1];
            vzR[2] = vR[2];
            vzL[1] = vL[1];
            vzL[2] = vL[2];
            bzR[1] = bR[1] * suRm;
            bzR[2] = bR[2] * suRm;
            bzL[1] = bL[1] * suLm;
            bzL[2] = bL[2] * suLm;
        }

        // Сгладим Скорость     Я добавил, не было изначально !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // vzR[1] = UZ[2] / UZ[0]
        // vzR[2] = UZ[3] / UZ[0]
        // vzL[1] = vzR[1]
        // vzL[2] = vzR[2]

        double UZL[8];
        double UZR[8];

        UZL[0] = rzL;
        UZL[4] = ezL;
        UZR[0] = rzR;
        UZR[4] = ezR;

        for (size_t ik = 0; ik < 3; ik++)
        {
            UZL[ik + 1] = vzL[ik] * rzL;
            UZL[ik + 5] = bzL[ik];
            UZR[ik + 1] = vzR[ik] * rzR;
            UZR[ik + 5] = bzR[ik];
        }

        Eigen::Vector3d qv;
        Eigen::Vector3d qb;

        if (SL >= wv)
        {
            qqq[0] = FL[0] - wv * UL[0];
            qqq[4] = FL[4] - wv * UL[4];

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] - wv * UL[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] - wv * UL[ik];
            }
        }
        else if (SL <= wv && SM >= wv)
        {
            qqq[0] = FL[0] + SL * (rzL - r1) - wv * UZL[0];
            qqq[4] = FL[4] + SL * (ezL - e1) - wv * UZL[4];

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FL[ik] + SL * (UZL[ik] - UL[ik]) - wv * UZL[ik];
            }
        }
        else if (SM <= wv && SR > wv)
        {
            qqq[0] = FR[0] + SR * (rzR - r2) - wv * UZR[0];
            qqq[4] = FR[4] + SR * (ezR - e2) - wv * UZR[4];

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] + SR * (UZR[ik] - UR[ik]) - wv * UZR[ik];
            }
        }
        else if (SR <= wv)
        {
            qqq[0] = FR[0] - wv * UR[0];
            qqq[4] = FR[4] - wv * UR[4];

            for (size_t ik = 1; ik < 4; ik++)
            {
                qv[ik - 1] = FR[ik] - wv * UR[ik];
            }

            for (size_t ik = 5; ik < 8; ik++)
            {
                qb[ik - 5] = FR[ik] - wv * UR[ik];
            }
        }
        else
        {
            cout << "error  3423167452  " << endl;
            for (short unsigned int ik = 0; ik < 8; ik++)
            {
                cout << "ik: " << qqq1[ik] << "   " << qqq2[ik] << ";" << endl;
            }
            whach(SR);
            whach(SL);
            whach(SM);
            whach(wv);
            exit(-1);
        }

        // Bn
        double SN = max(fabs(SL), fabs(SR));

        double wbn = 0.0;
        if (wv > SR)
        {
            wbn = wv * bR[0];
        }
        else if(wv < SL)
        {
            wbn = wv * bL[0];
        }
        else
        {
            wbn = wv * (bL[0] + bR[0]) / 2.0;
        }

        qb[0] = -SN * (bR[0] - bL[0]) - wbn;

        if (null_bn == true) qb[0] = 0.0;


        for (short unsigned int ik = 0; ik < 3; ik++)
        {
            qqq[ik + 1] = n(ik) * qv(0) + t(ik) * qv(1) + m(ik) * qv(2);
            qqq[ik + 5] = n(ik) * qb(0) + t(ik) * qb(1) + m(ik) * qb(2);
            qqq[ik + 5] = spi4 * qqq[ik + 5];
        }

        return;
    }

    return;
}

void Phys_param::raspad_testing(void)
{
    std::vector<double> qqq, qqq1, qqq2;
    qqq.resize(8);
    qqq1.resize(8);
    qqq2.resize(8);

    std::vector<double> konvect_left, konvect_right, konvect;
    PrintOptions Option = PrintOptions{};

    double dsr, dsc, dsl;

    qqq1[0] = 1.0;
    qqq1[1] = 1.0;
    qqq1[2] = 1.0;
    qqq1[3] = 1.0;
    qqq1[4] = 1.0;
    qqq1[5] = 5.0;
    qqq1[6] = 1.0;
    qqq1[7] = 1.0;

    qqq2[0] = 1.0;
    qqq2[1] = 1.0;
    qqq2[2] = 1.0;
    qqq2[3] = 1.0;
    qqq2[4] = 1.0;
    qqq2[5] = 1.0;
    qqq2[6] = 1.0;
    qqq2[7] = 1.0;

    double w = 0.0;

    this->chlld(1, 1.0, 0.0, 0.0,
        w, qqq1, qqq2, qqq, false, 3,
        konvect_left, konvect_right, konvect, dsr, dsc, dsl,
        Option);

    for (size_t ff = 0; ff < 8; ff++)
    {
        whach(ff);
        whach(qqq[ff]);
    }
}
