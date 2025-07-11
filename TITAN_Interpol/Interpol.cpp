#include "Interpol.h"

Interpol::Interpol(string name)
{
    cout << "Start: Interpol" << endl;
	std::ifstream in(name, std::ios::binary);
	if (!in) {
		cout << "Error 1329764965  Can not open file to reading: " + name << endl;
		exit(-1);
	}

	// Файл должен храниться в виде:
	// Сначала список строк с названиями переменных
	// Потом все ячейки в формате: три координаты и значения этих переменных

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

    in.close();

    // Делаем триангуляцию
    this->Delone = new Delaunay(this->points.begin(), this->points.end());


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

bool Interpol::Get_param(const double& x, const double& y, const double& z, std::unordered_map<string, double>& parameters, const Cell_handle& prev_cell, Cell_handle& next_cell)
{
    Point query(x, y, z);
    Cell_handle containing_cell = this->Delone->locate(query, prev_cell);

    next_cell = containing_cell;

    if (this->Delone->is_infinite(containing_cell)) 
    {
        return false;
    }

    // Получаем вершины тетраэдра 
    const  Point& p0 = containing_cell->vertex(0)->point();
    const  Point& p1 = containing_cell->vertex(1)->point();
    const  Point& p2 = containing_cell->vertex(2)->point();
    const  Point& p3 = containing_cell->vertex(3)->point();

    size_t i0 = containing_cell->vertex(0)->info();
    size_t i1 = containing_cell->vertex(1)->info();
    size_t i2 = containing_cell->vertex(2)->info();
    size_t i3 = containing_cell->vertex(3)->info();

    // Вычисляем барицентрические координаты 
    auto coords = barycentric_coordinates(query, Tetrahedron(p0, p1, p2, p3));



    for (const auto& nn : this->param_names)
    {
        parameters[nn] = coords[0] * this->Cells[i0]->parameters[nn] +
            coords[1] * this->Cells[i1]->parameters[nn] +
            coords[2] * this->Cells[i2]->parameters[nn] +
            coords[3] * this->Cells[i3]->parameters[nn];
    }

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
