#include "Interpol.h"

Interpol::Interpol(string name)
{
    // Степени r в интераоляции
    this->stepen["rho"] = 2.0;
    this->stepen["p"] = 2.0 * (5.0 / 3.0);
    this->stepen["divV"] = 1.0;
    this->stepen["Bx"] = 2.0;
    this->stepen["By"] = 2.0;
    this->stepen["Bz"] = 2.0;


    cout << "Start: Interpol" << endl;
	std::ifstream in(name, std::ios::binary);
	if (!in) {
		cout << "Error 1329764965  Can   not open file to reading: " + name << endl;
		exit(-1);
	}

	// Файл должен храниться в виде:
	// Сначала список строк с названиями переменных
	// Потом все ячейки в формате: три координаты и значения этих переменных

    in.read(reinterpret_cast<char*>(&this->razriv), sizeof(bool));

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

    if (this->razriv == true)
    {
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

        // Считываем третью зону
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, c);

            this->points_3.push_back({ {a, b, c}, i });

            for (const auto& i : this->param_names)
            {
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                A->parameters[i] = a;
            }

            this->Cells_3.push_back(A);
        }

        // Считываем четвёртую зону
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, c);

            this->points_4.push_back({ {a, b, c}, i });

            for (const auto& i : this->param_names)
            {
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                A->parameters[i] = a;
            }

            this->Cells_4.push_back(A);
        }

        // Считываем пятую зону
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, c);

            this->points_5.push_back({ {a, b, c}, i });

            for (const auto& i : this->param_names)
            {
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                A->parameters[i] = a;
            }

            this->Cells_5.push_back(A);
        }

        // Считываем шестую зону
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        for (size_t i = 0; i < size; ++i)
        {
            double a, b, c;
            in.read(reinterpret_cast<char*>(&a), sizeof(a));
            in.read(reinterpret_cast<char*>(&b), sizeof(b));
            in.read(reinterpret_cast<char*>(&c), sizeof(c));
            auto A = new Int_point(a, b, c);

            this->points_6.push_back({ {a, b, c}, i });

            for (const auto& i : this->param_names)
            {
                double a;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                A->parameters[i] = a;
            }

            this->Cells_6.push_back(A);
        }
    }

    int test_i = 0;
    in.read(reinterpret_cast<char*>(&test_i), sizeof(int));
    if (test_i != 121)
    {
        cout << "ERROR 1erjgiuheiuth8737y4tg08e4hoiufg4e" << endl;
        exit(-1);
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

            for (const auto& i : this->param_names)
            {
                if (i == "zone_geo") continue;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                in.read(reinterpret_cast<char*>(&b), sizeof(b));

                A->parameters[i + "_L"] = a;
                A->parameters[i + "_R"] = b;
            }

            this->Cells_TS.push_back(A);
        }
    }

    test_i = 0;
    in.read(reinterpret_cast<char*>(&test_i), sizeof(int));
    if (test_i != 122)
    {
        cout << "ERROR 2erjgiuheiuth8737y4tg08e4hoiufg4e" << endl;
        exit(-1);
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

        size_t size_x;
        size_t size_phi;
        in.read(reinterpret_cast<char*>(&size_x), sizeof(size_x));
        in.read(reinterpret_cast<char*>(&size_phi), sizeof(size_phi));
        this->Cells_HP_2.resize(boost::extents[size_x][size_phi]);

        for (size_t i = 0; i < size_phi; ++i)
        {
            for (size_t j = 0; j < size_x; ++j)
            {
                double a, b, c;
                in.read(reinterpret_cast<char*>(&a), sizeof(a));
                in.read(reinterpret_cast<char*>(&b), sizeof(b));
                in.read(reinterpret_cast<char*>(&c), sizeof(c));

                auto A = new Int_point(a, b, 0.0);
                A->parameters["r"] = c;

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

                this->Cells_HP_2[j][i] = A;
            }
        }
    }

    test_i = 0;
    in.read(reinterpret_cast<char*>(&test_i), sizeof(int));
    if (test_i != 123)
    {
        cout << "ERROR 3erjgiuheiuth8737y4tg08e4hoiufg4e" << endl;
        exit(-1);
    }

    // Считываем BS
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

            this->point_BS.push_back({ {a, b}, i });

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

            this->Cells_BS.push_back(A);
        }
    }

    test_i = 0;
    in.read(reinterpret_cast<char*>(&test_i), sizeof(int));
    if (test_i != 124)
    {
        cout << "ERROR 4erjgiuheiuth8737y4tg08e4hoiufg4e" << endl;
        exit(-1);
    }


    in.close();

    // Делаем триангуляцию
    this->Delone_1 = new Delaunay(this->points_1.begin(), this->points_1.end());

    if (this->razriv == true)
    {
        this->Delone_2 = new Delaunay(this->points_2.begin(), this->points_2.end());
        this->Delone_3 = new Delaunay(this->points_3.begin(), this->points_3.end());
        this->Delone_4 = new Delaunay(this->points_4.begin(), this->points_4.end());
        this->Delone_5 = new Delaunay(this->points_5.begin(), this->points_5.end());
        this->Delone_6 = new Delaunay(this->points_6.begin(), this->points_6.end());
    }

    this->Delone_TS = new Delaunay2(this->point_TS.begin(), this->point_TS.end());
    this->Delone_BS = new Delaunay2(this->point_BS.begin(), this->point_BS.end());
    this->Delone_HP_1 = new Delaunay2(this->point_HP_1.begin(), this->point_HP_1.end());

    cout << "END: Interpol" << endl;
}

Interpol::~Interpol()
{
    cout << "Start ~Interpol()" << endl;
    for (auto it = this->Cells_HP_2.data(); it != this->Cells_HP_2.data() + this->Cells_HP_2.num_elements(); ++it)
    {
        delete* it;
    }

    // Функция для очистки вектора указателей
    auto clearVector = [](auto& vec) {
        for (auto ptr : vec) {
            delete ptr;
        }
        vec.clear();
        };

    // Очищаем все векторы с Int_point*
    clearVector(Cells_1);
    clearVector(Cells_2);
    clearVector(Cells_3);
    clearVector(Cells_4);
    clearVector(Cells_5);
    clearVector(Cells_6);
    clearVector(Cells_TS);
    clearVector(Cells_BS);
    clearVector(Cells_HP_1);


    // Удаляем объекты Delaunay
    if (Delone_1 != nullptr) delete Delone_1;
    if (Delone_2 != nullptr) delete Delone_2;
    if (Delone_3 != nullptr) delete Delone_3;
    if (Delone_4 != nullptr) delete Delone_4;
    if (Delone_5 != nullptr) delete Delone_5;
    if (Delone_6 != nullptr) delete Delone_6;
    if (Delone_TS != nullptr) delete Delone_TS;
    if (Delone_BS != nullptr) delete Delone_BS;
    if (Delone_HP_1 != nullptr) delete Delone_HP_1;


    // Обнуляем указатели (опционально, но хорошая практика)
    Delone_1 = nullptr;
    Delone_2 = nullptr;
    Delone_3 = nullptr;
    Delone_4 = nullptr;
    Delone_5 = nullptr;
    Delone_6 = nullptr;
    Delone_TS = nullptr;
    Delone_BS = nullptr;
    Delone_HP_1 = nullptr;


    this->param_names.clear();
    cout << "End ~Interpol()" << endl;
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
    std::unordered_map<string, double>& parameters)
{
    short int this_zone;
    std::array<Cell_handle, 6> prev_cell;
    std::array<Cell_handle, 6> next_cell;
    for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();
    return this->Get_param(x, y, z, parameters, prev_cell, next_cell, this_zone);
}

bool Interpol::Get_param(const double& x, const double& y, const double& z,
    std::unordered_map<string, double>& parameters, short int& this_zone)
{
    std::array<Cell_handle, 6> prev_cell;
    std::array<Cell_handle, 6> next_cell;
    for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();
    return this->Get_param(x, y, z, parameters, prev_cell, next_cell, this_zone);
}

bool Interpol::Get_param(const double& x, const double& y, const double& z,
    std::unordered_map<string, double>& parameters, const std::array<Cell_handle, 6>& prev_cell, 
    std::array<Cell_handle, 6>& next_cell)
{
    short int this_zone;
    return this->Get_param(x, y, z, parameters, prev_cell, next_cell, this_zone);
}


bool Interpol::Get_param(const double& x, const double& y, const double& z, 
    std::unordered_map<string, double>& parameters, const std::array<Cell_handle, 6>& prev_cell, std::array<Cell_handle, 6>& next_cell,
    short int& this_zone)
{
    //cout << "A0" << endl;
    short int my_zone = 0;   // В какой зоне интерполяции находится точка
    this_zone = 0;
    std::unordered_map<string, double> param;
    double R_TS, R_HP, R_BS;
    bool aa1;

    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) 
    {
        std::cerr << "Error 87695476hturtguhge584   " << x << ", " << y << ", " << z << ")" << std::endl;
        exit(-1);
    }

    double radius = norm2(x, y, z);
    double radius2 = norm2(0.0, y, z);
    if (radius < 0.1) return false;

    // Далее определяем интерполяционную зону

    if (x < this->L6)
    {
        this_zone = 0;
        my_zone = 6;
        goto know_zone;
    }
    //cout << "A1" << endl; 
    aa1 = this->Get_TS(x, y, z, param);
    if (aa1 == false)
    {
        cout << "Error 9857686573jgh" << endl;
        exit(-1);
    }
    R_TS = param["r"];
    if (R_TS >= radius)
    {
        this_zone = 1;
        my_zone = 1;
        goto know_zone;
    }
    //cout << "A2" << endl;
    aa1 = this->Get_HP(x, y, z, param);
    if (aa1 == false)
    {
        cout << "Error hjtyj56y745ytgt" << endl;
        exit(-1);
    }
    R_HP = param["r"];

    if (x >= 0.0 && radius <= R_HP)
    {
        this_zone = 2;
        my_zone = 2;
        goto know_zone;
    }

    if (x < 0.0 && radius2 <= R_HP)
    {
        this_zone = 2;
        my_zone = 2;
        goto know_zone;
    }

    if (x <= 0.0)
    {
        this_zone = 3;
        my_zone = 5;
        goto know_zone;
    }
    //cout << "A3" << endl;
    aa1 = this->Get_BS(x, y, z, param);
    if (aa1 == false)
    {
        cout << "Error 9857686573jgh" << endl;
        exit(-1);
    }
    R_BS = param["r"];

    if (radius <= R_BS)
    {
        this_zone = 3;
        my_zone = 3;
        goto know_zone;
    }
    else
    {
        this_zone = 4;
        my_zone = 4;
        goto know_zone;
    }

    // Определили интерполяционную зону, теперь делаем интерполяцию
know_zone:
    //cout << "A4" << endl;
    Point query(x, y, z);
    Cell_handle containing_cell;
    std::vector <Int_point*>* CCC = nullptr;
    //cout << "A02" << endl;
    if (my_zone == 1 || this->razriv == false)
    {
        CCC = &this->Cells_1;

        this->mut_Delone_1.lock();
        containing_cell = this->Delone_1->locate(query, prev_cell[0]);
        this->mut_Delone_1.unlock();

        if (this->Delone_1->is_infinite(containing_cell))
        {
            //cout << "Delone Infinit  " << x << " " << y << " " << z << endl;
            return false;
            // Находим ближайшую вершину
            this->mut_Delone_1.lock();
            Vertex_handle nearest_vertex = this->Delone_1->nearest_vertex(query);
            containing_cell = nearest_vertex->cell();  // Берём любой смежный треугольник
            query = nearest_vertex->point();
            this->mut_Delone_1.unlock();
        }
        next_cell[0] = containing_cell;
    }
    else if((my_zone == 2))
    {
        CCC = &this->Cells_2;
        containing_cell = this->Delone_2->locate(query, prev_cell[1]);
        if (this->Delone_2->is_infinite(containing_cell))
        {
            return false;
        }
        next_cell[1] = containing_cell;
    }
    else if ((my_zone == 3))
    {
        CCC = &this->Cells_3;
        containing_cell = this->Delone_3->locate(query, prev_cell[2]);
        if (this->Delone_3->is_infinite(containing_cell))
        {
            return false;
        }
        next_cell[2] = containing_cell;
    }
    else if ((my_zone == 4))
    {
        CCC = &this->Cells_4;
        containing_cell = this->Delone_4->locate(query, prev_cell[3]);
        if (this->Delone_4->is_infinite(containing_cell))
        {
            return false;
        }
        next_cell[3] = containing_cell;
    }
    else if ((my_zone == 5))
    {
        CCC = &this->Cells_5;
        containing_cell = this->Delone_5->locate(query, prev_cell[4]);
        if (this->Delone_5->is_infinite(containing_cell))
        {
            return false;
        }
        next_cell[4] = containing_cell;
    }
    else if ((my_zone == 6))
    {
        CCC = &this->Cells_6;
        containing_cell = this->Delone_6->locate(query, prev_cell[5]);
        if (this->Delone_6->is_infinite(containing_cell))
        {
            return false;
        }
        next_cell[5] = containing_cell;
    }
    //cout << "A5" << endl;
   

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
   //cout << "A6" << endl;
    // Определяем зону в общих областях интерполяции
    if (my_zone == 6 || my_zone == 5)
    {
        // Определяем текущую зону
        short int zone[5];
        for (short int i = 0; i < 5; i++) zone[i] = 0;

        short int kk = static_cast<short int>(std::round((*CCC)[i0]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        if (kk == 1) kk = 2;
        zone[kk]++;
        kk = static_cast<short int>(std::round((*CCC)[i1]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        if (kk == 1) kk = 2;
        zone[kk]++;
        kk = static_cast<short int>(std::round((*CCC)[i2]->parameters["zone_geo"]));
        if (kk <= 0) kk = 0;
        if (kk == 1) kk = 2;
        zone[kk]++;
        kk = static_cast<short int>(std::round((*CCC)[i3]->parameters["zone_geo"]));
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
    //cout << "A7" << endl;
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

   // cout << "A8" << endl;
    for (const auto& nn : this->param_names)
    {
        parameters[nn] = 0.0;
        short int kl = 0;
        for (const auto& ii : i_n)
        {
            parameters[nn] += coords[kl] * (*CCC)[ii]->parameters[nn] * pow(r[kl], this->stepen[nn]);
            kl++;
        }

        parameters[nn] /= pow(r[4], this->stepen[nn]);
    }
    
    //cout << "A9" << endl;
    return true;
}


bool Interpol::Get_real_cells(const double& x, const double& y, const double& z,
    vector<int>& num_cell, vector<double>& koeff_cell, const Cell_handle& prev_cell, Cell_handle& next_cell)
{
    if (this->razriv == true)
    {
        cout << "Error 9j4fuheuohfguweorghuiyesrhgesgr" << endl;
        exit(-1);
    }

    double radius = norm2(x, y, z);
    double radius2 = norm2(0.0, y, z);
    if (radius < 0.1) return false;

    // Определили интерполяционную зону, теперь делаем интерполяцию
    // 
    //cout << "A01" << endl;
    Point query(x, y, z);
    Cell_handle containing_cell;
    std::vector <Int_point*>* CCC = nullptr;
    //cout << "A02" << endl;
    
    CCC = &this->Cells_1;

    this->mut_Delone_1.lock();
    containing_cell = this->Delone_1->locate(query, prev_cell);
    this->mut_Delone_1.unlock();

    if (this->Delone_1->is_infinite(containing_cell))
    {
        return false;
    }
    next_cell = containing_cell;
    
    //cout << "A03" << endl;
    // Получаем вершины тетраэдра 
    Point& p0 = containing_cell->vertex(0)->point();
    Point& p1 = containing_cell->vertex(1)->point();
    Point& p2 = containing_cell->vertex(2)->point();
    Point& p3 = containing_cell->vertex(3)->point();

    size_t i0 = containing_cell->vertex(0)->info();
    size_t i1 = containing_cell->vertex(1)->info();
    size_t i2 = containing_cell->vertex(2)->info();
    size_t i3 = containing_cell->vertex(3)->info();

    //cout << "A031" << endl;

    if (num_cell.size() < 4) num_cell.resize(4);
    if (koeff_cell.size() < 4) koeff_cell.resize(4);

    num_cell[0] = i0;
    num_cell[1] = i1;
    num_cell[2] = i2;
    num_cell[3] = i3;

    //cout << "A04" << endl;

    // Вычисляем барицентрические координаты 
    auto coords = barycentric_coordinates(query, Tetrahedron(p0, p1, p2, p3));

    //cout << "A05" << endl;

    koeff_cell[0] = coords[0];
	koeff_cell[1] = coords[1];
	koeff_cell[2] = coords[2];
	koeff_cell[3] = coords[3];

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


double bilinearInterpolationInRectangle(
    double& x1, double& y1, double& f1,  // левый нижний
    double& x2, double& y2, double& f2,  // правый нижний  
    double& x3, double& y3, double& f3,  // правый верхний
    double& x4, double& y4, double& f4,  // левый верхний
    double& x0, double& y0) {           // целевая точка

    // Для прямоугольника: x1 = x4, x2 = x3, y1 = y2, y3 = y4
    // Проверяем что это действительно прямоугольник
    if (x1 != x4 || x2 != x3 || y1 != y2 || y3 != y4) {
        std::cerr << "Warning: Not a perfect rectangle!" << std::endl;
        exit(-1);
    }

    // Нормированные координаты внутри прямоугольника
    double u = (x0 - x1) / (x2 - x1);
    double v = (y0 - y1) / (y4 - y1);

    // Ограничиваем параметры [0, 1] если точка за границами
    u = std::max(0.0, std::min(1.0, u));
    v = std::max(0.0, std::min(1.0, v));

    // Билинейная интерполяция
    double bottom = f1 * (1 - u) + f2 * u;    // интерполяция по нижней стороне
    double top = f4 * (1 - u) + f3 * u;       // интерполяция по верхней стороне

    return bottom * (1 - v) + top * v;        // интерполяция по вертикали
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
        if (x < this->L6)
        {
            cout << "HP net pri x < " << this->L6 << endl;
            return false;
        }

        //cout << "C1" << endl;
        double r_1, x_1, phi_1;

        r_1 = sqrt(y * y + z * z);
        x_1 = x;
        phi_1 = polar_angle(y, z);

        if (r_1 < 0.001) phi_1 = 0.0;

        Point2 query(x_1, phi_1);
        double x0 = x_1;
        double y0 = phi_1;

        short int i1, j1;
        short int i2, j2;

        // Количество строк (первое измерение)
        size_t rows = this->Cells_HP_2.shape()[0];

        // Количество столбцов (второе измерение)  
        size_t cols = this->Cells_HP_2.shape()[1];

        i1 = rows - 1;
        j1 = cols - 1;

        //cout << "C2 " << rows << " " << cols << " " << i1 << " " << j1 << endl;

        while (this->Cells_HP_2[i1][j1]->center[0] > x0)
        {
            i1--;
            if (i1 < 0)
            {
                i1 = 0;
                break;
            }
        }

        //cout << "C3 " << rows << " " << cols << " " << i1 << " " << j1 << endl;

        while (this->Cells_HP_2[i1][j1]->center[1] > y0)
        {
            j1--;
            if (j1 < 0)
            {
                j1 = 0;
                break;
            }
        }

        //cout << "C4 " << rows << " " << cols << " " << i1 << " " << j1 << endl;

        i2 = i1 + 1;
        j2 = j1 + 1;
        if (i2 >= rows) i2 = rows - 1;
        if (j2 >= cols) j2 = cols - 1;

        //cout << "C5 " << rows << " " << cols << " " << i1 << " " << j1 << endl;

        double x1 = this->Cells_HP_2[i1][j1]->center[0];
        double x2 = this->Cells_HP_2[i2][j1]->center[0];
        double x3 = x2;
        double x4 = x1;

        double y1 = this->Cells_HP_2[i1][j1]->center[1];
        double y2 = y1;
        double y3 = this->Cells_HP_2[i1][j2]->center[1];
        double y4 = y3;

        for (auto& [key, _] : this->Cells_HP_2[i1][j1]->parameters)
        {
            double f1 = this->Cells_HP_2[i1][j1]->parameters[key];
            double f2 = this->Cells_HP_2[i2][j1]->parameters[key];
            double f3 = this->Cells_HP_2[i2][j2]->parameters[key];
            double f4 = this->Cells_HP_2[i1][j2]->parameters[key];

            parameters[key] = bilinearInterpolationInRectangle(x1, y1, f1, x2, y2, f2,
                x3, y3, f3, x4, y4, f4, x0, y0);
        }


        return true;
    }
}

bool Interpol::Get_BS(const double& x, const double& y, const double& z,
    std::unordered_map<string, double>& parameters)
{
    if (x < 0.0)
    {
        cout << "BS net pri x < " << 0.0 << endl;
        return false;
    }

    if (true)
    {
        double x_ = x;
        if (x_ < 0.0001) x_ = 0.0001;

        double r_1, the_1, phi_1;

        r_1 = sqrt(x_ * x_ + y * y + z * z);
        the_1 = acos(x_ / r_1);
        phi_1 = polar_angle(y, z);

        Point2 query(the_1, phi_1);

        //cout << "S1 " << endl;

        Face_handle face = this->Delone_BS->locate(query);
        if (this->Delone_BS->is_infinite(face))
        {
            // Находим ближайшую вершину
            Vertex_handle2 nearest_vertex = this->Delone_BS->nearest_vertex(query);
            face = nearest_vertex->face();  // Берём любой смежный треугольник
            Point2 p0 = face->vertex(0)->point();
            size_t idx0 = face->vertex(0)->info();

            for (auto& [key, _] : this->Cells_BS[idx0]->parameters)
            {
                double f0 = this->Cells_BS[idx0]->parameters[key];
                parameters[key] = f0;
            }
            return true;
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
        for (auto& [key, _] : this->Cells_BS[idx0]->parameters)
        {
            double f0 = this->Cells_BS[idx0]->parameters[key];
            double f1 = this->Cells_BS[idx1]->parameters[key];
            double f2 = this->Cells_BS[idx2]->parameters[key];

            // Преобразуем CGAL::FT в double
            parameters[key] = CGAL::to_double(alpha) * f0 +
                CGAL::to_double(beta) * f1 +
                CGAL::to_double(gamma) * f2;
        }
        //cout << "S4" << endl;
        //cout << "F" << endl;
        return true;
    }
}