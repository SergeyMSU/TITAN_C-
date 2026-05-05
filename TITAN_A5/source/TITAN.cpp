// TITAN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "Header.h"
using namespace std;
//class Setka;

// Структура для хранения точки в полярных координатах
struct PolarPoint 
{
    double phi; // угол в радианах [0, pi]
    double r;   // радиус
};


// Функция чтения данных из файла и преобразования в полярные координаты
std::vector<PolarPoint> readAndConvert(const std::string& filename) 
{
    std::vector<PolarPoint> points;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error open file: " << filename << std::endl;
        return points;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Пропуск пустых строк
        if (line.empty()) continue;

        std::istringstream iss(line);
        double x, y;
        if (!(iss >> x >> y)) 
        {
            std::cerr << "Error read: " << line << std::endl;
            continue;
        }

        // Преобразование в полярные координаты
        double r = std::sqrt(x * x + y * y);
        double phi = polar_angle(x, y);

        points.push_back({ phi, r });
    }

    // Проверка, что углы идут по возрастанию (для надёжности)
    for (size_t i = 1; i < points.size(); ++i) 
    {
        if (points[i].phi < points[i - 1].phi) {
            std::cerr << "Ugli ne monotonni " << i << std::endl;
        }
    }

    return points;
}

// Функция линейной интерполяции r по заданному углу phi
// Используется последовательный перебор для поиска интервала (т.к. данных немного)
double interpolateR(const std::vector<PolarPoint>& points, double phi_query) {
    
    if (points.empty()) 
    {
        std::cerr << "Ошибка: массив точек пуст." << std::endl;
        return 0.0;
    }


    // Если запрос меньше первого угла — возвращаем r первой точки
    if (phi_query <= points.front().phi) {
        return points.front().r;
    }

    // Если запрос больше или равен последнему углу — возвращаем r последней точки
    if (phi_query >= points.back().phi) {
        return points.back().r;
    }

    // Последовательный поиск интервала, содержащего phi_query
    for (size_t i = 0; i < points.size() - 1; ++i)
    {
        double phi1 = points[i].phi;
        double phi2 = points[i + 1].phi;

        // Нашли интервал [phi1, phi2], в который попадает phi_query
        if (phi_query >= phi1 && phi_query <= phi2) {
            double r1 = points[i].r;
            double r2 = points[i + 1].r;

            // Линейная интерполяция: r = r1 + (r2 - r1) * (phi_query - phi1) / (phi2 - phi1)
            double t = (phi_query - phi1) / (phi2 - phi1);
            return r1 + t * (r2 - r1);
        }
    }

    // Сюда не должны попасть при корректных данных, но на всякий случай:
    std::cerr << "Ошибка: не удалось найти интервал для phi = " << phi_query << std::endl;
    return 0.0;
}


int main()
{
    cout << "Start Programm" << endl;




    if (false)
    {
        Eigen::Vector3d vec, cc, vv;
        Eigen::Vector3d vec2, cc2, vv2;
        Eigen::Matrix3d Matr;              // Матрица перехода из HGI в мои
        Eigen::Matrix3d Matr2;             // Матрица перехода моих в HGI
        Eigen::Matrix3d Matr3;             // Матрица перехода из Эклиптических координат в мои

        Matr << -0.9958639688067077, 0.01776569097515556, 0.08910295088675518,
            0.07561695085992419, 0.7057402284561812, 0.7044237408557894,
            -0.05036896241933166, 0.7082479157489926, -0.7041646522383864;

        Matr2 << -0.9958639688067080, 0.0756169508599243, -0.0503689624193315,
            0.0177656909751554, 0.7057402284561816, 0.7082479157489927,
            0.0891029508867553, 0.7044237408557898, -0.7041646522383865;



        vec << -0.9952, -0.0669, -0.0709;
        vv << 1.0, 0.0, 0.0;

        cc = Matr * vec;

        cout << cc[0] << " " << cc[1] << " " << cc[2] << endl;
        double cosAngle = cc.dot(vv) / (cc.norm() * vv.norm());
        cosAngle = std::clamp(cosAngle, -1.0, 1.0);
        cout << std::acos(cosAngle) * 180.0 / const_pi << endl;

        Matr3 << -0.2510319412434562, -0.9637264652805977, 0.0906325801977802,
            -0.5738471722585876, 0.2235689862406300, 0.7878555269096995,
            -0.7795398561756577, 0.1457676524785322, -0.6091546635498517;


        vec2 << -0.173236, -0.982923, -0.062060;
        cc = Matr3 * vec2;
        cout << cc[0] << " " << cc[1] << " " << cc[2] << endl;
        cosAngle = cc.dot(vv) / (cc.norm() * vv.norm());
        cosAngle = std::clamp(cosAngle, -1.0, 1.0);
        cout << std::acos(cosAngle) * 180.0 / const_pi << endl;

        return 0;
    }

    // Создаём основную сетку из файлов вспомогательных сеток
    //Setka S1 = Setka("SDK1_2D_Setka.bin", "SDK1_krug_setka.bin", 60);
    Setka S1 = Setka("SDK_A5.2_2D_Setka.bin", "SDK1_krug_setka.bin", 60);
    //Setka S1 = Setka("SDK2_2D_Setka.bin", "SDK1_krug_setka.bin", 60);

    S1.geo->L6 = -150.0; // Этот параметр важен при загрузке сетки так как определяет до куда выделять HP
    S1.geo->L7 = -300.0;
    S1.geo->tetta1 = 2.89;
    S1.geo->tetta2 = 2.64;



    // Обязательный блок настройки основной сетки
    if (true)
    {
        // Считаем старый файл поверхностей что-бы приблизительно подвинуть их в нужное место
        S1.Read_old_surface("ASurf_Save00591.bin");

        // Теперь передвигаем сетку к поверхностям
        S1.Move_to_surf(S1.Surf1);

        // Автоматически подстраиваем геометрические параметры сетки (сгущение и т.д.) под новые поверхности
        S1.auto_set_luch_geo_parameter(0);

        // Считаем объёмы, площади и другие геометрические характеристики
        S1.Calculating_measure(0);
        S1.Calculating_measure(1);

        // Задаём граничные грани
        S1.Init_boundary_grans();
    }



    // Если надо интерполируем значения из другой сетки
    //S1.PereInterpolate("For_intertpolate_0081-no_razriv-with_MK.bin", true);

    // Считываем физические параметры и геометрическое положение узлов из файла (предыдущего расчёта)
    //S1.Download_cell_parameters("parameters_0060.bin"); 
    //  

    //S1.Download_cell_parameters("parameters_0079.bin"); 
    //   c 10 начал ручное передвижение сетки, потом в 11 его подкорректировал, с 12 начал считать


    //S1.Download_cell_parameters("parameters_promeg_119.bin");  //60   
    //S1.Download_cell_parameters("parameters_A5_0076.bin");  //60   
    //S1.Download_cell_parameters("parameters_A5_0085.bin");  //60   
    // на 46 подвинул поверхности вручную
    // 
    // 63 - Лакс и вторым порядком + особый снос в сверхзвуке
    // 60 - запустил движение поверхностей
    // 68 - без движения, везде HLLD
    // 74 - поменял L6 на -190 (было -200)
    // 
    //S1.Download_cell_parameters("parameters_A5_0047.bin");  
    S1.Download_cell_parameters("parameters_A5-2_0001.bin");  
    

    //S1.Download_cell_parameters("parameters_promeg_1121.bin");  


    //S1.Download_cell_parameters("parameters_0219.bin");    
    // S1.Download_cell_parameters("parameters_promeg_1112.bin");
    //S1.Download_cell_parameters("parameters_promeg_1118.bin");
    //S1.Download_cell_parameters("parameters_promeg_1112.bin");

    //S1.Download_cell_parameters("parameters_0219.bin");


    S1.geo->L6 = -150.0;   // Этот параметр важен при загрузке сетки
    S1.geo->L7 = -300.0;
    S1.geo->tetta1 = 2.89;
    S1.geo->tetta2 = 2.64;

    //S1.PereInterpolate("For_intertpolate_0082-no_razriv-with_MK.bin", true);


    //S1.Download_cell_parameters("parameters_promeg_1124.bin");

    // Ещё один блок обязательной настройки
    if (true)
    {
        // Точно задаём положение внутренней границы сетки
        S1.geo->R0 = S1.phys_param->R_0;
        S1.geo->R1 = 20.0; // 20.0;

        // Автоматически подстраиваем геометрические параметры сетки под новые положения узлов
        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Luch)
        {
            i->dvigenie(0);
        }

        // Инициализируем TVD (находим соседей и т.д.)
        S1.Init_TVD();

        S1.Find_Yzel_Sosed_for_sglag();

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }
    }


    // Ручное движение TS
    if (false)
    {
        cout << "Hand TS move" << endl;
        // 1. Считываем файл и преобразуем в полярные координаты
        std::string filename = "TS_.txt"; // укажите правильный путь к файлу
        std::vector<PolarPoint> polarPoints = readAndConvert(filename);

        cout << "TS 1: " << interpolateR(polarPoints, 0.0) << endl;
        cout << "TS 2: " << interpolateR(polarPoints, 1.0) << endl;
        cout << "TS 3: " << interpolateR(polarPoints, 2.0) << endl;
        cout << "TS 4: " << interpolateR(polarPoints, 3.0) << endl;

        for (auto& i : S1.All_Yzel)
        {
            if (i->type != Type_yzel::TS) continue;

            double x = i->coord[0][0];
            double y = i->coord[0][1];
            double z = i->coord[0][2];
            double r = norm2(x, y, z);
            double phi = polar_angle(x, norm2(0.0, y, z));
            double r_interp = interpolateR(polarPoints, phi);

            i->coord[1][0] *= (r + 1.0 * (r_interp - r)) / r;
            i->coord[1][1] *= (r + 1.0 * (r_interp - r)) / r;
            i->coord[1][2] *= (r + 1.0 * (r_interp - r)) / r;

            i->coord[0][0] = i->coord[1][0];
            i->coord[0][1] = i->coord[1][1];
            i->coord[0][2] = i->coord[1][2];
        }

        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        cout << "END Hand TS move" << endl;
    }

    // Ручное движение HP
    if (false)
    {
        cout << "Hand HP move" << endl;
        // 1. Считываем файл и преобразуем в полярные координаты
        std::string filename = "HP_.txt"; // укажите правильный путь к файлу
        std::vector<PolarPoint> polarPoints = readAndConvert(filename);


        for (auto& i : S1.All_Yzel)
        {
            if (i->type != Type_yzel::HP) continue;
            double x = i->coord[0][0];
            double y = i->coord[0][1];
            double z = i->coord[0][2];
            double r = norm2(x, y, z);
            double phi = polar_angle(x, norm2(0.0, y, z));
            double r_interp = interpolateR(polarPoints, phi);


            i->coord[1][0] *= (r + 1.0 * (r_interp - r)) / r;
            i->coord[1][1] *= (r + 1.0 * (r_interp - r)) / r;
            i->coord[1][2] *= (r + 1.0 * (r_interp - r)) / r;

            i->coord[0][0] = i->coord[1][0];
            i->coord[0][1] = i->coord[1][1];
            i->coord[0][2] = i->coord[1][2];
        }

        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        cout << "END Hand HP move" << endl;

        // Двигаем узлы на невыделяемой части HP
        int now2 = 0;
        short int NN = S1.D_Luch[0].size() - 1;
        int kk = 10;
        for (auto& L : S1.D_Luch)
        {
            /*int kk = 1;
            while (true)
            {
                kk++;
                if (L[S1.geo->N4 - kk]->Yzels_opor[0]->coord[now2][0] > S1.geo->L6) break;
            }*/

            double h1 = norm2(0.0, L[S1.geo->N4 - kk]->Yzels_opor[1]->coord[now2][1],
                L[S1.geo->N4 - kk]->Yzels_opor[1]->coord[now2][2]);

            cout << "h1 = " << h1 << endl;


            for (short int i = S1.geo->N4 - kk + 1; i <= NN; i++)
            {
                double h = h1;
                auto yz = L[i]->Yzels_opor[1];
                double hh = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
                yz->coord[now2][1] = yz->coord[now2][1] * h / hh;
                yz->coord[now2][2] = yz->coord[now2][2] * h / hh;
            }
        }

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
    }

    //  Ручное изменение BS
    if (false)
    {
        cout << "Hand BS move" << endl;
        // 1. Считываем файл и преобразуем в полярные координаты
        std::string filename = "BS_.txt"; // укажите правильный путь к файлу
        std::vector<PolarPoint> polarPoints = readAndConvert(filename);


        for (auto& i : S1.All_Yzel)
        {
            if (i->type != Type_yzel::BS) continue;

            //cout << "x do = " << i->coord[0][0] << endl;
            double x = i->coord[0][0];
            double y = i->coord[0][1];
            double z = i->coord[0][2];
            double r = norm2(x, y, z);
            double phi = polar_angle(x, norm2(0.0, y, z));
            double r_interp = interpolateR(polarPoints, phi);

            //cout << "DO phi = " << phi * 180.0/const_pi << "   r do = " << r << "     r posle = " << r_interp << endl;

            i->coord[1][0] *= (r_interp / r);
            i->coord[1][1] *= (r_interp / r);
            i->coord[1][2] *= (r_interp / r);

            i->coord[0][0] = i->coord[1][0];
            i->coord[0][1] = i->coord[1][1];
            i->coord[0][2] = i->coord[1][2];

            //cout << "x posle = " << i->coord[0][0] << endl;
        }


        S1.auto_set_luch_geo_parameter(0);


        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        cout << "END Hand BS move" << endl;

        // Надо подвинуть узлы, которые продолжают BS

        short int NN = S1.A_Luch[0].size() - 1;
        int now2 = 0;
        for (short int i = 0; i < S1.A_Luch.size(); i++)
        {
            auto yz = S1.A_Luch[i][NN]->Yzels_opor[3];
            double H = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
            for (auto& L : S1.B_Luch[i])
            {
                auto yyz = L->Yzels_opor[3];
                double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
                yyz->coord[now2][1] *= H / h2;
                yyz->coord[now2][2] *= H / h2;
            }
            for (auto& L : S1.E_Luch[i])
            {
                auto yyz = L->Yzels_opor[2];
                double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
                yyz->coord[now2][1] *= H / h2;
                yyz->coord[now2][2] *= H / h2;
            }
            for (auto& L : S1.D_Luch[i])
            {
                auto yyz = L->Yzels_opor[2];
                double h2 = norm2(0.0, yyz->coord[now2][1], yyz->coord[now2][2]);
                yyz->coord[now2][1] *= H / h2;
                yyz->coord[now2][2] *= H / h2;
            }
        }

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
    }

    // Ручное движение HP
    if (false)
    {
        cout << "Hand HP move" << endl;
        // 1. Считываем файл и преобразуем в полярные координаты
        std::string filename = "HP_.txt"; // укажите правильный путь к файлу
        std::vector<PolarPoint> polarPoints = readAndConvert(filename);


        for (auto& i : S1.All_Yzel)
        {
            if (i->type != Type_yzel::HP) continue;
            double x = i->coord[0][0];
            double y = i->coord[0][1];
            double z = i->coord[0][2];
            double r = norm2(x, y, z);
            double phi = polar_angle(x, norm2(0.0, y, z));
            double r_interp = interpolateR(polarPoints, phi);


            i->coord[1][0] *= (r + 1.0 * (r_interp - r)) / r;
            i->coord[1][1] *= (r + 1.0 * (r_interp - r)) / r;
            i->coord[1][2] *= (r + 1.0 * (r_interp - r)) / r;

            i->coord[0][0] = i->coord[1][0];
            i->coord[0][1] = i->coord[1][1];
            i->coord[0][2] = i->coord[1][2];
        }

        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        cout << "END Hand HP move" << endl;

        // Двигаем узлы на невыделяемой части HP
        int now2 = 0;
        short int NN = S1.D_Luch[0].size() - 1;
        int kk = 10;
        for (auto& L : S1.D_Luch)
        {
            /*int kk = 1;
            while (true)
            {
                kk++;
                if (L[S1.geo->N4 - kk]->Yzels_opor[0]->coord[now2][0] > S1.geo->L6) break;
            }*/

            double h1 = norm2(0.0, L[S1.geo->N4 - kk]->Yzels_opor[1]->coord[now2][1],
                L[S1.geo->N4 - kk]->Yzels_opor[1]->coord[now2][2]);

            for (short int i = S1.geo->N4 - kk + 1; i <= NN; i++)
            {
                double h = h1;
                auto yz = L[i]->Yzels_opor[1];
                double hh = norm2(0.0, yz->coord[now2][1], yz->coord[now2][2]);
                yz->coord[now2][1] = yz->coord[now2][1] * h / hh;
                yz->coord[now2][2] = yz->coord[now2][2] * h / hh;
            }
        }

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.auto_set_luch_geo_parameter(0);

        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }

        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
    }


 

    
    //S1.PereInterpolate("For_intertpolate_0084-no_razriv-with_MK.bin", false, false);

    //S1.Smooth_head_TS();


    cout << "A-" << endl;
    // Задаём начальные и граничные условия
    //S1.Init_physics();
    

    // Блок начальной визуализации сетки для проверки корректного построения
    if (true)
    {
        S1.Tecplot_print_cell_plane_parameters();
        S1.Tecplot_print_all_lush_in_2D();
        S1.Tecplot_print_2D_setka(0.0, 0.0, 1.0, -0.00001, "init_setka_2d_(0, 0, 1, 0)_");
        //S1.Tecplot_print_2D_setka(0.0, 1.0, 0.0, -0.00001, "init_setka_2d_(0, 1, 0, 0)_");
        //S1.Tecplot_print_2D_setka(0.0, 1.0, 1.0, -0.00001, "init_setka_2d_(0, 1, 1, 0)_");
        S1.Tecplot_print_all_gran_in_surface("TS");
        S1.Tecplot_print_all_gran_in_surface("HP");
        S1.Tecplot_print_all_gran_in_surface("BS");
        S1.Tecplot_print_plane_lush(30);
    }

    // Выбор основного алгоритма расчёта (в данной функции представлены все варианты расчёта: атомы, мгд и т.д.), см. саму функцию
    //S1.Algoritm(23, &S1);

    
    //S1.Algoritm(101, &S1);
    //S1.Algoritm(23, &S1);
    //S1.Algoritm(10, &S1);
    //S1.Algoritm(23, &S1);
    //S1.Algoritm(101, &S1);

    //S1.Algoritm(10, &S1);
    //S1.Algoritm(5, &S1);ts

    S1.cooling = new CoolingFunction();
    S1.heating = new CoolingFunction();

    S1.cooling->ReadCoolingFunction("combined_cooling_function.txt");
    S1.heating->ReadCoolingFunction("combined_heating_function.txt");

    
    //S1.Algoritm(1, &S1);

    //S1.Write_file_for_FCMHD();
    S1.Read_file_for_FCMHD();

    //Dust_spectra DDD = Dust_spectra();



    //return 0;

    /// Далее следует всё, что касается визуализации сетки



    if (false)
    {
        // Планировал запустить дальше перестройку сорта 2, потом зоны 2, 4, 6
        //S1.Algoritm(2);
        //S1.Algoritm(8);
        //S1.Algoritm(5);
        //S1.Print_fH(4, Type_Gran_surf::BS, 1.0, 0.0, 0.0, 5.0 * const_pi/180.0);
        //S1.Print_fH(4, Type_Gran_surf::HP, 1.0, 0.0, 0.0, 5.0 * const_pi / 180.0);
        //S1.Print_fH(2, Type_Gran_surf::TS, 1.0, 0.0, 0.0, 5.0 * const_pi / 180.0);

        //S1.Print_f_proect_in_gran(1);
        //S1.Print_f_proect_in_gran(2);
        //S1.Print_f_proect_in_gran(3);

        /*S1.Print_f_proect_in_cell(10.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(15.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(20.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(25.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(30.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(35.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(40.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(70.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(80.0, 0.0, 0.0);
        S1.Print_f_proect_in_cell(45.0, 0.0, 0.0);*/

        cout << "AABB" << endl;
        /*S1.Print_SpSm(17.0, 0.0, 0.0);
        S1.Print_SpSm(20.0, 0.0, 0.0);
        S1.Print_SpSm(25.0, 0.0, 0.0);
        S1.Print_SpSm(1.0, 0.0, 0.0);
        S1.Print_SpSm(5.0, 0.0, 0.0);
        S1.Print_SpSm(10.0, 0.0, 0.0);
        S1.Print_SpSm(15.0, 0.0, 0.0);*/

        /*S1.Print_pui(17.0, 0.0, 0.0);
        S1.Print_pui(20.0, 0.0, 0.0);
        S1.Print_pui(25.0, 0.0, 0.0);
        S1.Print_pui(1.0, 0.0, 0.0);
        S1.Print_pui(5.0, 0.0, 0.0);
        S1.Print_pui(10.0, 0.0, 0.0);
        S1.Print_pui(15.0, 0.0, 0.0);
        S1.Print_pui(28.0, 0.0, 0.0);
        S1.Print_pui(50.0, 0.0, 0.0);
        S1.Print_pui(100.0, 0.0, 0.0);
        S1.Print_pui(200.0, 0.0, 0.0);*/

        return 0;
    }


    //S1.Save_cell_parameters("parameters_A5-2_0001.bin");
    //S1.Save_cell_parameters("parameters_0138.bin");
    //S1.Save_cell_pui_parameters("parameters_0026.bin");

    /*S1.Edges_create();
    S1.Culc_divergence_in_cell();
    S1.Culc_gradient_in_cell();
    S1.Culc_rotors_in_cell();
    S1.Culc_usual_rotors_in_cell();*/

    if (false)
    {
        // Надо улучшить ротеры вблизи разрывов
        for (auto& gr : S1.Gran_TS)
        {
            auto C1 = gr->cells[1];
            auto C2 = gr->cells_TVD[1];

            C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
            C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
            C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];

            C1->parameters[0]["gradBB_x"] = C2->parameters[0]["gradBB_x"];
            C1->parameters[0]["gradBB_y"] = C2->parameters[0]["gradBB_y"];
            C1->parameters[0]["gradBB_z"] = C2->parameters[0]["gradBB_z"];

            C1 = gr->cells[0];
            C2 = gr->cells_TVD[0];

            C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
            C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
            C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];

            C1->parameters[0]["gradBB_x"] = C2->parameters[0]["gradBB_x"];
            C1->parameters[0]["gradBB_y"] = C2->parameters[0]["gradBB_y"];
            C1->parameters[0]["gradBB_z"] = C2->parameters[0]["gradBB_z"];
        }

        for (auto& gr : S1.Gran_HP)
        {
            auto C1 = gr->cells[0];
            auto C2 = gr->cells_TVD[0];

            C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
            C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
            C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];

            C1->parameters[0]["gradBB_x"] = C2->parameters[0]["gradBB_x"];
            C1->parameters[0]["gradBB_y"] = C2->parameters[0]["gradBB_y"];
            C1->parameters[0]["gradBB_z"] = C2->parameters[0]["gradBB_z"];


            C1 = gr->cells[0];
            C2 = gr->cells_TVD[0];

            C1->parameters[0]["rotB_x"] = C2->parameters[0]["rotB_x"];
            C1->parameters[0]["rotB_y"] = C2->parameters[0]["rotB_y"];
            C1->parameters[0]["rotB_z"] = C2->parameters[0]["rotB_z"];

            C1->parameters[0]["gradBB_x"] = C2->parameters[0]["gradBB_x"];
            C1->parameters[0]["gradBB_y"] = C2->parameters[0]["gradBB_y"];
            C1->parameters[0]["gradBB_z"] = C2->parameters[0]["gradBB_z"];
        }


    }

    S1.Save_for_interpolate("For_intertpolate_0085-no_razriv-with_MK.bin", false);
    //return 0;
    Interpol SS = Interpol("For_intertpolate_0085-no_razriv-with_MK.bin");

    //S1.Save_for_interpolate("For_intertpolate_0059-.bin", false);
    //Interpol SS = Interpol("For_intertpolate_0059-.bin");

    cout << "AAA" << endl;

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0), "_(1, 0, 0)_", 500.0);

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(-1.0, 0.0, 0.0), "_(-1, 0, 0)_", 500.0);

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0), "_(0, 1, 0)_", 500.0);

    S1.Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");
    S1.Tecplot_print_2D(&SS, 0.0, 1.0, 0.0, -0.00001, "_2d_(0, 1, 0, 0)_");


    Eigen::Vector3d eex(1.0, 0.0, 0.0);
    Eigen::Vector3d eey(0.0, -sqrt(2.0) / 2.0, sqrt(2.0) / 2.0);
    Eigen::Vector3d centr_sys(0.0, 0.0, 0.0);
    S1.Tecplot_print_2D(&SS, 0.0, 1.0, 1.0, -0.00001, "_2d_(0, 1, 1, 0)_", false, eex, eey, centr_sys);


    cout << "F " << endl;

    S1.Tecplot_print_cell_plane_parameters();


    cout << "YSPEX" << endl;


    //S1.Tecplot_print_all_yzel_in_3D("SDK1");
    
    

    //S1.Tecplot_print_krug_yzel_in_3D(1);
    //S1.Tecplot_print_krug_yzel_in_3D(2);

    //S1.Tecplot_print_all_lush_in_2D();
    //S1.Tecplot_print_All_surfase_in_2D();
    //S1.Tecplot_print_plane_lush(0);
    //S1.Tecplot_print_plane_surfase(0);
    //S1.Tecplot_print_all_gran_in_cell();
    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");
    //S1.Tecplot_print_all_yzel_with_condition();
}

