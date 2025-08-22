#include "Interpol.h"

Interpol::Interpol(string name)
{
    // Степени r в интераоляции
    this->stepen["rho"] = 2.0;
    this->stepen["p"] = 2.0 * (5.0 / 3.0);
    this->stepen["divV"] = 1.0;
    this->stepen["Bx"] = 1.0;
    this->stepen["By"] = 1.0;
    this->stepen["Bz"] = 1.0;


    cout << "Start: Interpol" << endl;
	std::ifstream in(name, std::ios::binary);
	if (!in) {
		cout << "Error 1329764965  Can not open file to reading: " + name << endl;
		exit(-1);
	}

	// Файл должен храниться в виде:
	// Сначала список строк с названиями переменных
	// Потом все ячейки в формате: три координаты и значения этих переменных

    in.read(reinterpret_cast<char*>(&this->L6), sizeof(this->L6));

     // Читаем количество строк
    size_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(size));

    // Читаем каждую строку
    for (size_t i = 0; i < size; ++i) 
    {
        // Читаем длину строки
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        // Читаем саму строку
        std::string str(str_size, '\0');
        in.read(&str[0], str_size);
        this->param_names.push_back(str);
    }

    cout << "All parameters: " << endl;
    for (const auto& i : this->param_names)
    {
        cout << i << "  ";
    }
    cout << endl;

    // Считываем первую зону
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    for (size_t i = 0; i < size; ++i)
    {
        double a, b, c;
        in.read(reinterpret_cast<char*>(&a), sizeof(a));
        in.read(reinterpret_cast<char*>(&b), sizeof(b));
        in.read(reinterpret_cast<char*>(&c), sizeof(c));
        auto A = new Int_point(a, b, c);

        this->points_1.push_back({ {a, b, c}, i });

        for (const auto& i : this->param_names)
        {
            double a;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            A->parameters[i] = a;
        }

        this->Cells_1.push_back(A);
    }

    // Считываем вторую зону
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    for (size_t i = 0; i < size; ++i)
    {
        double a, b, c;
        in.read(reinterpret_cast<char*>(&a), sizeof(a));
        in.read(reinterpret_cast<char*>(&b), sizeof(b));
        in.read(reinterpret_cast<char*>(&c), sizeof(c));
        auto A = new Int_point(a, b, c);

        this->points_2.push_back({ {a, b, c}, i });

        for (const auto& i : this->param_names)
        {
            double a;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            A->parameters[i] = a;
        }

        this->Cells_2.push_back(A);
    }


    // Считываем TS
    if (true)
    {
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, 0.0);
            A->parameters["r"] = c;

            this->point_TS.push_back({ {a, b}, i });

            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            A->parameters["nx"] = a;
            A->parameters["ny"] = b;
            A->parameters["nz"] = c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));

            A->parameters["rho_L"] = a;
            A->parameters["rho_R"] = b;

            this->Cells_TS.push_back(A);
        }
    }

    // Считываем HP
    if (true)
    {
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, 0.0);
            A->parameters["r"] = c;

            this->point_HP_1.push_back({ {a, b}, i });

            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            A->parameters["nx"] = a;
            A->parameters["ny"] = b;
            A->parameters["nz"] = c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));

            A->parameters["rho_L"] = a;
            A->parameters["rho_R"] = b;

            this->Cells_HP_1.push_back(A);
        }

        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, 0.0);
            A->parameters["r"] = c;

            this->point_HP_2.push_back({ {a, b}, i });

            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            A->parameters["nx"] = a;
            A->parameters["ny"] = b;
            A->parameters["nz"] = c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));

            A->parameters["rho_L"] = a;
            A->parameters["rho_R"] = b;

            this->Cells_HP_2.push_back(A);
        }
    }




    in.close();

    // Делаем триангуляцию
    this->Delone_1 = new Delaunay(this->points_1.begin(), this->points_1.end());
    this->Delone_2 = new Delaunay(this->points_2.begin(), this->points_2.end());

    this->Delone_TS = new Delaunay2(this->point_TS.begin(), this->point_TS.end());
    this->Delone_HP_1 = new Delaunay2(this->point_HP_1.begin(), this->point_HP_1.end());

    cout << "END: Interpol" << endl;
}

Interpol::~Interpol()
{
    delete Delone_1;
    for (auto& i : this->Cells_1)
    {
        delete i;
    }

    this->Cells_1.clear();
    this->points_1.clear();
    this->param_names.clear();
}

// Вычисление объёма тетраэдра
double tetrahedron_volume(const Point& A, const Point& B, const Point& C, const Point& D) {
    Vector AB = B - A;
    Vector AC = C - A;
    Vector AD = D - A;
    return CGAL::scalar_product(CGAL::cross_product(AB, AC), AD) / 6.0;
}

std::array<double, 4> barycentric_coordinates(const Point& P, const Tetrahedron& tet) {
    const Point& A = tet.vertex(0);
    const Point& B = tet.vertex(1);
    const Point& C = tet.vertex(2);
    const Point& D = tet.vertex(3);

    double V = tetrahedron_volume(A, B, C, D);
    double lambda1 = tetrahedron_volume(P, B, C, D) / V;
    double lambda2 = tetrahedron_volume(A, P, C, D) / V;
    double lambda3 = tetrahedron_volume(A, B, P, D) / V;
    double lambda4 = 1.0 - lambda1 - lambda2 - lambda3;

    return { lambda1, lambda2, lambda3, lambda4};
}

bool Interpol::Get_param(const double& x, const double& y, const double& z,
    std::unordered_map<string, double>& parameters, const Cell_handle& prev_cell, Cell_handle& next_cell)
{
    short int this_zone;
    return this->Get_param(x, y, z, parameters, prev_cell, next_cell, this_zone);
}

bool Interpol::Get_param(const double& x, const double& y, const double& z, 
    std::unordered_map<string, double>& parameters, const Cell_handle& prev_cell, Cell_handle& next_cell, 
    short int& this_zone)
{
    //cout << "A0" << endl;

    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        std::cerr << "Ошибка: недопустимые координаты (" << x << ", " << y << ", " << z << ")" << std::endl;
        exit(-1);
    }
    double radius = norm2(x, y, z);

    if (radius < 0.1) return false;

    std::unordered_map<string, double> param;
    bool aa1 = this->Get_TS(x, y, z, param);
    if (aa1 == false)
    {
        cout << "Error 9857686573jgh" << endl;
        exit(-1);
    }

    //cout << "A01" << endl;

    this_zone = 0;
    if (param["r"] >= radius)
    {
        this_zone = 1;
    }

    Point query(x, y, z);
    Cell_handle containing_cell;
    std::vector <Int_point*>* CCC;
    //cout << "A02" << endl;
    if (this_zone == 1)
    {
        CCC = &this->Cells_1;
        containing_cell = this->Delone_1->locate(query);
        if (this->Delone_1->is_infinite(containing_cell))
        {
            return false;
        }
    }
    else
    {
        CCC = &this->Cells_2;
        containing_cell = this->Delone_2->locate(query);
        if (this->Delone_2->is_infinite(containing_cell))
        {
            return false;
        }
    }
    //cout << "A1" << endl;
    next_cell = containing_cell;

    //cout << "A2" << endl;
    // Получаем вершины тетраэдра 
    Point& p0 = containing_cell->vertex(0)->point();
    Point& p1 = containing_cell->vertex(1)->point();
    Point& p2 = containing_cell->vertex(2)->point();
    Point& p3 = containing_cell->vertex(3)->point();

    size_t i0 = containing_cell->vertex(0)->info();
    size_t i1 = containing_cell->vertex(1)->info();
    size_t i2 = containing_cell->vertex(2)->info();
    size_t i3 = containing_cell->vertex(3)->info();

    // Вычисляем барицентрические координаты 
    auto coords = barycentric_coordinates(query, Tetrahedron(p0, p1, p2, p3));

    //cout << "A2" << endl;

    if (this_zone == 0)
    {
        // Определяем текущую зону
        short int zone[5];
        if (true)
        {
            for (short int i = 0; i < 5; i++) zone[i] = 0;

            short int kk = static_cast<short int>(std::round(this->Cells_2[i0]->parameters["zone_geo"]));
            if (kk <= 0) kk = 0;
            if (kk == 1) kk = 2;
            zone[kk]++;
            kk = static_cast<short int>(std::round(this->Cells_2[i1]->parameters["zone_geo"]));
            if (kk <= 0) kk = 0;
            if (kk == 1) kk = 2;
            zone[kk]++;
            kk = static_cast<short int>(std::round(this->Cells_2[i2]->parameters["zone_geo"]));
            if (kk <= 0) kk = 0;
            if (kk == 1) kk = 2;
            zone[kk]++;
            kk = static_cast<short int>(std::round(this->Cells_2[i3]->parameters["zone_geo"]));
            if (kk <= 0) kk = 0;
            if (kk == 1) kk = 2;
            zone[kk]++;

            kk = 0;
            for (short int i = 1; i < 5; i++)
            {
                if (zone[i] > kk)
                {
                    this_zone = i;
                    kk = zone[i];
                }
            }

            if (this_zone <= 1)
            {
                cout << "Error 9867540975" << endl;
                exit(-1);
            }

            if (this_zone == 2 || this_zone == 3)
            {
                double Q = coords[0] * this->Cells_2[i0]->parameters["Q"] +
                    coords[1] * this->Cells_2[i1]->parameters["Q"] +
                    coords[2] * this->Cells_2[i2]->parameters["Q"] +
                    coords[3] * this->Cells_2[i3]->parameters["Q"];
                double rho = coords[0] * this->Cells_2[i0]->parameters["rho"] +
                    coords[1] * this->Cells_2[i1]->parameters["rho"] +
                    coords[2] * this->Cells_2[i2]->parameters["rho"] +
                    coords[3] * this->Cells_2[i3]->parameters["rho"];
                if (Q / rho < 70.0)
                {
                    this_zone = 2;
                }
                else
                {
                    this_zone = 3;
                }
            }


        }
    }

    //cout << "A3" << endl;
    vector<size_t> i_n(4);
    i_n[0] = i0;
    i_n[1] = i1;
    i_n[2] = i2;
    i_n[3] = i3;

    double r[5];
    if (this_zone == 1)
    {
        r[0] = norm2(p0[0], p0[1], p0[2]);
        r[1] = norm2(p1[0], p1[1], p1[2]);
        r[2] = norm2(p2[0], p2[1], p2[2]);
        r[3] = norm2(p3[0], p3[1], p3[2]);
        r[4] = norm2(x, y, z);
        if (r[4] < 0.01)
        {
            r[0] = 1.0;
            r[1] = 1.0;
            r[2] = 1.0;
            r[3] = 1.0;
            r[4] = 1.0;
        }
    }
    else
    {
        r[0] = 1.0;
        r[1] = 1.0;
        r[2] = 1.0;
        r[3] = 1.0;
        r[4] = 1.0;
    }

    //cout << "A4" << endl;
    for (const auto& nn : this->param_names)
    {
        /*parameters[nn] = coords[0] * this->Cells[i0]->parameters[nn] +
            coords[1] * this->Cells[i1]->parameters[nn] +
            coords[2] * this->Cells[i2]->parameters[nn] +
            coords[3] * this->Cells[i3]->parameters[nn];*/

        parameters[nn] = 0.0;
        short int kl = 0;
        for (const auto& ii : i_n)
        {
            parameters[nn] += coords[kl] * (*CCC)[ii]->parameters[nn] * pow(r[kl], this->stepen[nn]);
            kl++;
        }

        parameters[nn] /= pow(r[4], this->stepen[nn]);
    }
    
    //cout << "A5" << endl;
    return true;
}

// Функция для вычисления барицентрических координат
void compute_barycentric(const Point2& a, const Point2& b, const Point2& c,
    const Point2& p, FT& alpha, FT& beta, FT& gamma)
{
    FT denom = (b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y());
    alpha = ((b.y() - c.y()) * (p.x() - c.x()) + (c.x() - b.x()) * (p.y() - c.y())) / denom;
    beta = ((c.y() - a.y()) * (p.x() - c.x()) + (a.x() - c.x()) * (p.y() - c.y())) / denom;
    gamma = FT(1) - alpha - beta;
}


bool Interpol::Get_TS(const double& x, const double& y, const double& z,
    std::unordered_map<string, double>& parameters)
{
    //cout << "S1 " << endl;
    double r_1, the_1, phi_1;

    r_1 = sqrt(x * x + y * y + z * z);
    the_1 = acos(x / r_1);
    phi_1 = polar_angle(y, z);

    //cout << "S2  " << the_1 << "    "
     //    << phi_1 << endl;

    Point2 query(the_1, phi_1);

    //cout << "S3 " << endl;

    Face_handle face = this->Delone_TS->locate(query);
    if (this->Delone_TS->is_infinite(face)) 
    {
        // Находим ближайшую вершину
        Vertex_handle2 nearest_vertex = this->Delone_TS->nearest_vertex(query);
        face = nearest_vertex->face();  // Берём любой смежный треугольник
    }
    //cout << "B" << endl;
    // Получаем вершины треугольника
    Point2 p0 = face->vertex(0)->point();
    Point2 p1 = face->vertex(1)->point();
    Point2 p2 = face->vertex(2)->point();
    //cout << "C" << endl;
    FT alpha, beta, gamma;
    compute_barycentric(p0, p1, p2, query, alpha, beta, gamma);
    //cout << "D" << endl;
    size_t idx0 = face->vertex(0)->info(); // Номер точки p0
    size_t idx1 = face->vertex(1)->info(); // Номер точки p1
    size_t idx2 = face->vertex(2)->info(); // Номер точки p2

    // Интерполяция параметров
    for (auto& [key, _] : this->Cells_TS[idx0]->parameters) 
    {
        double f0 = this->Cells_TS[idx0]->parameters[key];
        double f1 = this->Cells_TS[idx1]->parameters[key];
        double f2 = this->Cells_TS[idx2]->parameters[key];

        // Преобразуем CGAL::FT в double
        parameters[key] = CGAL::to_double(alpha) * f0 +
            CGAL::to_double(beta) * f1 +
            CGAL::to_double(gamma) * f2;
    }
    //cout << "F" << endl;
    return true;
}


// Функция для решения системы линейных уравнений 2x2
bool solveLinearSystem2x2(double a11, double a12, double a21, double a22,
    double b1, double b2, double& x, double& y) {
    double det = a11 * a22 - a12 * a21;

    if (std::abs(det) < 1e-12) {
        return false; // Вырожденная система
    }

    x = (b1 * a22 - b2 * a12) / det;
    y = (a11 * b2 - a21 * b1) / det;

    return true;
}

// Билинейная интерполяция в произвольном четырёхугольнике
double bilinearInterpolationInQuadrilateral(
    const Point2& p1, double f1,  // (x1, y1), f1
    const Point2& p2, double f2,  // (x2, y2), f2  
    const Point2& p3, double f3,  // (x3, y3), f3
    const Point2& p4, double f4,  // (x4, y4), f4
    const Point2& p0) {           // (x0, y0) - целевая точка

    // Преобразуем в декартовы координаты
    double x1 = CGAL::to_double(p1.x()), y1 = CGAL::to_double(p1.y());
    double x2 = CGAL::to_double(p2.x()), y2 = CGAL::to_double(p2.y());
    double x3 = CGAL::to_double(p3.x()), y3 = CGAL::to_double(p3.y());
    double x4 = CGAL::to_double(p4.x()), y4 = CGAL::to_double(p4.y());
    double x0 = CGAL::to_double(p0.x()), y0 = CGAL::to_double(p0.y());

    // Параметрическое представление четырёхугольника
    // p(u,v) = (1-u)(1-v)*p1 + u(1-v)*p2 + uv*p3 + (1-u)v*p4

    // Находим параметры u, v такие что:
    // x0 = (1-u)(1-v)*x1 + u(1-v)*x2 + u*v*x3 + (1-u)*v*x4
    // y0 = (1-u)(1-v)*y1 + u(1-v)*y2 + u*v*y3 + (1-u)*v*y4

    // Перепишем в виде системы уравнений:
    // A*u + B*v + C*u*v = D
    // E*u + F*v + G*u*v = H

    double A = x2 - x1;
    double B = x4 - x1;
    double C = x1 - x2 - x4 + x3;
    double D = x0 - x1;

    double E = y2 - y1;
    double F = y4 - y1;
    double G = y1 - y2 - y4 + y3;
    double H = y0 - y1;

    // Решаем систему методом Ньютона
    double u = 0.5, v = 0.5; // Начальное приближение
    const int max_iterations = 50;
    const double tolerance = 1e-8;

    for (int iter = 0; iter < max_iterations; ++iter) {
        double fu = A * u + B * v + C * u * v - D;
        double fv = E * u + F * v + G * u * v - H;

        if (std::abs(fu) < tolerance && std::abs(fv) < tolerance) {
            break;
        }

        // Матрица Якоби
        double J11 = A + C * v;
        double J12 = B + C * u;
        double J21 = E + G * v;
        double J22 = F + G * u;

        double detJ = J11 * J22 - J12 * J21;

        if (std::abs(detJ) < 1e-12) {
            // Вырожденная матрица Якоби, используем запасной метод
            u = (x0 - x1) / (x2 - x1); // Простая линейная аппроксимация
            v = (y0 - y1) / (y4 - y1);
            break;
        }

        // Коррекция Ньютона
        double du = (-fu * J22 + fv * J12) / detJ;
        double dv = (fu * J21 - fv * J11) / detJ;

        u += du;
        v += dv;

        // Ограничиваем параметры [0,1]
        u = std::max(0.0, std::min(1.0, u));
        v = std::max(0.0, std::min(1.0, v));
    }

    // Интерполируем значение
    double result = (1 - u) * (1 - v) * f1 + u * (1 - v) * f2 + u * v * f3 + (1 - u) * v * f4;

    return result;
}


bool Interpol::Get_HP(const double& x, const double& y, const double& z,
    std::unordered_map<string, double>& parameters)
{
    // ОБРАТИТЕ ВНИМАНИЕ, что для x > 0 и для x < 0 используются различные виды интерполяции
    if (x >= 0.0)
    {
        double r_1, the_1, phi_1;

        r_1 = sqrt(x * x + y * y + z * z);
        the_1 = acos(x / r_1);
        phi_1 = polar_angle(y, z);

        Point2 query(the_1, phi_1);

        //cout << "S1 " << endl;

        Face_handle face = this->Delone_HP_1->locate(query);
        if (this->Delone_HP_1->is_infinite(face))
        {
            // Находим ближайшую вершину
            Vertex_handle2 nearest_vertex = this->Delone_HP_1->nearest_vertex(query);
            face = nearest_vertex->face();  // Берём любой смежный треугольник
        }
        //cout << "S2" << endl;
        // Получаем вершины треугольника
        Point2 p0 = face->vertex(0)->point();
        Point2 p1 = face->vertex(1)->point();
        Point2 p2 = face->vertex(2)->point();
        //cout << "C" << endl;
        FT alpha, beta, gamma;
        compute_barycentric(p0, p1, p2, query, alpha, beta, gamma);
        //cout << "D" << endl;
        size_t idx0 = face->vertex(0)->info(); // Номер точки p0
        size_t idx1 = face->vertex(1)->info(); // Номер точки p1
        size_t idx2 = face->vertex(2)->info(); // Номер точки p2
        //cout << "S3" << endl;
        // Интерполяция параметров
        for (auto& [key, _] : this->Cells_HP_1[idx0]->parameters)
        {
            double f0 = this->Cells_HP_1[idx0]->parameters[key];
            double f1 = this->Cells_HP_1[idx1]->parameters[key];
            double f2 = this->Cells_HP_1[idx2]->parameters[key];

            // Преобразуем CGAL::FT в double
            parameters[key] = CGAL::to_double(alpha) * f0 +
                CGAL::to_double(beta) * f1 +
                CGAL::to_double(gamma) * f2;
        }
        //cout << "S4" << endl;
        //cout << "F" << endl;
        return true;
    }
    else
    {
        cout << "SS00" << endl;
        if (x < this->L6)
        {
            cout << "HP net pri x < " << this->L6 << endl;
            return false;
        }

        double r_1, x_1, phi_1;

        r_1 = sqrt(y * y + z * z);
        x_1 = x;
        phi_1 = polar_angle(y, z);

        Point2 query(x_1, phi_1);
        double x0 = x_1;
        double y0 = phi_1;

        size_t l00, l10, l11, l01;
        double x00 = -10000.0, y00 = -100000.0;
        double x10 = 10000.0, y10 = -100000.0;
        double x11 = 10000.0, y11 = 100000.0;
        double x01 = -10000.0, y01 = 100000.0;
        Point2 p00, p10, p11, p01;
        // Инициализируем переменные для поиска ближайших точек в каждом квадранте
        
        cout << "SS0" << endl;
        for (const auto& point_pair : this->point_HP_2) 
        {
            Point2 p = point_pair.first;

            // 1. Левая нижняя точка: x <= x0, y <= y0
            if (p.x() <= x0 && p.y() <= y0 && y00 > y0)
            {
                x00 = norm2(0.0, x0 - p.x(), y0 - p.y());
                l00 = point_pair.second;
                p00 = p;
            }


            // 2. Правая нижняя точка: x >= x0, y <= y0
            if (p.x() >= x0 && p.y() <= y0 && norm2(0.0, x0 - p.x(), y0 - p.y()) < x10)
            {
                x10 = norm2(0.0, x0 - p.x(), y0 - p.y());
                l10 = point_pair.second;
                p10 = p;
            }

            // 3. Правая верхняя точка: x >= x0, y >= y0
            if (p.x() >= x0 && p.y() >= y0 && norm2(0.0, x0 - p.x(), y0 - p.y()) < x11)
            {
                x11 = norm2(0.0, x0 - p.x(), y0 - p.y());
                l11 = point_pair.second;
                p11 = p;
            }

            // 4. Левая верхняя точка: x <= x0, y >= y0
            if (p.x() <= x0 && p.y() >= y0 && norm2(0.0, x0 - p.x(), y0 - p.y()) < x01)
            {
                x01 = norm2(0.0, x0 - p.x(), y0 - p.y());
                l01 = point_pair.second;
                p01 = p;
            }
        }

        cout << p00.x() << " " << p00.y() << endl;
        cout << p10.x() << " " << p10.y() << endl;
        cout << p11.x() << " " << p11.y() << endl;
        cout << p01.x() << " " << p01.y() << endl;
        cout << query.x() << " " << query.y() << endl << endl;

        cout << "SS1" << endl;
        for (auto& [key, _] : this->Cells_HP_1[l00]->parameters)
        {
            double f00 = this->Cells_HP_1[l00]->parameters[key];
            double f10 = this->Cells_HP_1[l10]->parameters[key];
            double f11 = this->Cells_HP_1[l11]->parameters[key];
            double f01 = this->Cells_HP_1[l01]->parameters[key];

            parameters[key] = bilinearInterpolationInQuadrilateral(p00, f00, p10, f10,
                p11, f11, p01, f01, query);
            //cout << f00 << " " << f10 << " " << f11 << " " << f01 << "  = " << parameters[key] << endl;
        }
        cout << "SS2" << endl;
        return true;
    }
}