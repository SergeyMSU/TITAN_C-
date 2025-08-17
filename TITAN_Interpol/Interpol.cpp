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

    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    for (size_t i = 0; i < size; ++i)
    {
        double a, b, c;
        in.read(reinterpret_cast<char*>(&a), sizeof(a));
        in.read(reinterpret_cast<char*>(&b), sizeof(b));
        in.read(reinterpret_cast<char*>(&c), sizeof(c));
        auto A = new Int_point(a, b, c);

        this->points.push_back({ {a, b, c}, i });

        for (const auto& i : this->param_names)
        {
            double a;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            A->parameters[i] = a;
        }

        this->Cells.push_back(A);
    }


    // Считываем центральную ячейку
    if (true)
    {
        double a, b, c;
        in.read(reinterpret_cast<char*>(&a), sizeof(a));
        in.read(reinterpret_cast<char*>(&b), sizeof(b));
        in.read(reinterpret_cast<char*>(&c), sizeof(c));
        auto A = new Int_point(a, b, c);

        this->points.push_back({ {a, b, c}, size});

        for (const auto& i : this->param_names)
        {
            double a;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            A->parameters[i] = a;
        }

        this->Cells.push_back(A);
    }

    // Считываем TS
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


    // Считываем точки TS- и TS+
    if (true)
    {
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, c);

            this->points_TS_1.push_back({ {a, b, c}, i });

            for (const auto& i : this->param_names)
            {
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                A->parameters[i] = a;
            }

            this->Cells_TS_1.push_back(A);
        }

        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, c);

            this->points_TS_2.push_back({ {a, b, c}, i });

            for (const auto& i : this->param_names)
            {
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                A->parameters[i] = a;
            }

            this->Cells_TS_2.push_back(A);
        }

    }



    in.close();

    // Делаем триангуляцию
    this->Delone = new Delaunay(this->points.begin(), this->points.end());
    this->Delone_TS_1 = new Delaunay(this->points_TS_1.begin(), this->points_TS_1.end());
    this->Delone_TS_2 = new Delaunay(this->points_TS_2.begin(), this->points_TS_2.end());

    this->Delone_TS = new Delaunay2(this->point_TS.begin(), this->point_TS.end());

    cout << "END: Interpol" << endl;
}

Interpol::~Interpol()
{
    delete Delone;
    for (auto& i : this->Cells)
    {
        delete i;
    }

    this->Cells.clear();
    this->points.clear();
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
    Point query(x, y, z);
    Cell_handle containing_cell = this->Delone->locate(query, prev_cell);
    //cout << "A1" << endl;
    next_cell = containing_cell;

    if (this->Delone->is_infinite(containing_cell)) 
    {
        return false;
    }
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
    //cout << "A3" << endl;
    // Вычисляем барицентрические координаты 
    auto coords = barycentric_coordinates(query, Tetrahedron(p0, p1, p2, p3));
    //cout << "A4" << endl;

    // Определяем текущую зону
    this_zone = 0;
    short int zone[5];
    if (true)
    {
        for (short int i = 0; i < 5; i++) zone[i] = 0;

        short int kk = static_cast<short int>(std::round(this->Cells[i0]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        zone[kk]++;
        kk = static_cast<short int>(std::round(this->Cells[i1]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        zone[kk]++;
        kk = static_cast<short int>(std::round(this->Cells[i2]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        zone[kk]++;
        kk = static_cast<short int>(std::round(this->Cells[i3]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        zone[kk]++;
    }

    std::vector <Int_point*>* CCC;
    CCC = &this->Cells;
    bool in_TS = false;

    //cout << "A5" << endl;
    // В этом случае мы вблизи TS
    if (zone[1] > 0 && zone[2] > 0)
    {
        Cell_handle containing_cell2;
        in_TS = true;
        std::unordered_map<string, double> param;
        bool aa1 = this->Get_TS(x, y, z, param);
        if (aa1 == false)
        {
            cout << "Error 9857686573jgh" << endl;
            exit(-1);
        }
        double radius = norm2(x, y, z);
        if (radius < param["r"])
        {
            this_zone = 1;
            containing_cell2 = this->Delone_TS_1->locate(query);
            if (this->Delone_TS_1->is_infinite(containing_cell2))
            {
                cout << "Error ewty456y36746" << endl;
                exit(-1);
            }
            CCC = &this->Cells_TS_1;
        }
        else
        {
            this_zone = 2;
            containing_cell2 = this->Delone_TS_2->locate(query);
            if (this->Delone_TS_2->is_infinite(containing_cell2))
            {
                cout << "Error r6u568u57854yhgft" << endl;
                exit(-1);
            }
            CCC = &this->Cells_TS_2;
        }
        p0 = containing_cell2->vertex(0)->point();
        p1 = containing_cell2->vertex(1)->point();
        p2 = containing_cell2->vertex(2)->point();
        p3 = containing_cell2->vertex(3)->point();

        i0 = containing_cell2->vertex(0)->info();
        i1 = containing_cell2->vertex(1)->info();
        i2 = containing_cell2->vertex(2)->info();
        i3 = containing_cell2->vertex(3)->info();

        // Вычисляем барицентрические координаты 
        coords = barycentric_coordinates(query, Tetrahedron(p0, p1, p2, p3));
    }
    else
    {
        short int kk = 0;
        for (short int i = 1; i < 5; i++)
        {
            if (zone[i] > kk)
            {
                this_zone = i;
                kk = zone[i];
            }
        }

        if (this_zone == 2 || this_zone == 3)
        {
            double Q = coords[0] * (*CCC)[i0]->parameters["Q"] +
                coords[1] * (*CCC)[i1]->parameters["Q"] +
                coords[2] * (*CCC)[i2]->parameters["Q"] +
                coords[3] * (*CCC)[i3]->parameters["Q"];
            double rho = coords[0] * (*CCC)[i0]->parameters["rho"] +
                coords[1] * (*CCC)[i1]->parameters["rho"] +
                coords[2] * (*CCC)[i2]->parameters["rho"] +
                coords[3] * (*CCC)[i3]->parameters["rho"];
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
    //cout << "B" << endl;

    //cout << "A6" << endl;
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

    //cout << "A7" << endl;
    //cout << (*CCC).size() << " " << i0 << " " << i1 << " " << i2 << " " << i3 << endl;
    //cout << in_TS << endl;
    //cout << this->Cells.size() << endl;
    //cout << this->Cells_TS_1.size() << endl;
    //cout << this->Cells_TS_2.size() << endl;

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
    
    //cout << "A8" << endl;
    /*cout << query[0] << " " << query[1] << " " << query[2] << endl;
    cout << p0[0] << " " << p0[1] << " " << p0[2] << endl;
    cout << p1[0] << " " << p1[1] << " " << p1[2] << endl;
    cout << p2[0] << " " << p2[1] << " " << p2[2] << endl;
    cout << p3[0] << " " << p3[1] << " " << p3[2] << endl;
    cout << parameters["rho"] << endl;
    cout << this->Cells[i0]->parameters["rho"] << endl;
    cout << this->Cells[i1]->parameters["rho"] << endl;
    cout << this->Cells[i2]->parameters["rho"] << endl;
    cout << this->Cells[i3]->parameters["rho"] << endl;

    exit(-1);*/

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
    //cout << "A" << endl;
    double r_1, the_1, phi_1;

    r_1 = sqrt(x * x + y * y + z * z);
    the_1 = acos(x / r_1);
    phi_1 = polar_angle(y, z);

    Point2 query(the_1, phi_1);
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