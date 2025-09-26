// TITAN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "Header.h"
using namespace std;
//class Setka;

int main()
{
    cout << "Start Programm" << endl;
    
    //auto phys_param = new Phys_param();
    //phys_param->raspad_testing();
    //return 0;

    Setka S1 = Setka();

    S1.Read_old_surface("ASurf_Save00591.bin");
    S1.Move_to_surf(S1.Surf1);

    S1.auto_set_luch_geo_parameter(0);
    cout << "A " << endl;
    S1.Calculating_measure(0);
    cout << "B " << endl;
    S1.Calculating_measure(1);
    cout << "B2 " << endl;
    
    S1.Init_boundary_grans();
    cout << "C " << endl;

    //S1.Download_cell_parameters("parameters_0137.bin");   // 107   119
    //S1.Download_cell_parameters("parameters_0057.bin");   // 107
    S1.Download_cell_parameters("parameters_0208.bin");   // 107

    // 19 стартовая точка от которой две параллели с пикапами и без
    // 32 с пикапами
    // 62 включи TVD
    // 76 перед тем, как отключить все ТВД и все костыли

    S1.geo->R0 = S1.phys_param->R_0;

    // 23 полностью установленное решение без Пикапов (у контакта есть артефакт нужно сглаживание по
    // углу увеличить)

    cout << "C2 " << endl;

    S1.auto_set_luch_geo_parameter(0);

    //S1.Smooth_head_HP2(); // Ручное сглаживание
    //S1.Smooth_HP1(); // Ручное сглаживание

    S1.Init_TVD();
    cout << "D2 " << endl;

    S1.Init_physics();

    cout << "E " << endl;

    //S1.Smooth_head_TS2();
    //S1.Smooth_head_HP2();

    S1.Tecplot_print_cell_plane_parameters();
    S1.Tecplot_print_all_lush_in_2D();
    S1.Tecplot_print_all_cell_in_3D();

    //S1.Print_SpSm(20.0, 0.0, 0.0);
    //S1.Print_SpSm(10.0, 0.0, 0.0);
    //S1.Print_SpSm(40.0, 0.0, 0.0);
    //return 0;

    if (false)
    {
        S1.Algoritm(5);
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
        S1.Print_pui(200.0, 0.0, 0.0);

        return 0;*/
    }


    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");

    

    S1.Find_Yzel_Sosed_for_BS();

    S1.Smooth_angle_HP();
    S1.Smooth_head_HP3();
    S1.Smooth_head_TS3();




    for (int i = 1; i <= 6 * 22; i++) // 6 * 2
    {
        auto start = std::chrono::high_resolution_clock::now();
        cout << "IIIII = " << i << endl;

        S1.Go(false, 400, 1); // 400   1
        S1.Go(true, 100, 1); // 400   1 
        S1.Smooth_head_HP3();
        S1.Smooth_head_TS3();

        S1.Tecplot_print_cell_plane_parameters();
        S1.Tecplot_print_all_lush_in_2D();
        S1.Tecplot_print_all_gran_in_surface("TS");
        S1.Tecplot_print_all_gran_in_surface("HP");
        S1.Tecplot_print_all_gran_in_surface("BS");
        //S1.Go(true, 100, 1);
        //S1.Tecplot_print_cell_plane_parameters();

        //S1.Init_physics();

        if (i % 12 == 0)
        {
            string namn = "parameters_promeg_11" + to_string(i) + ".bin";
            S1.Save_cell_parameters(namn);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << "Execution time: " << duration.count()/1000.0/60.0 << " minutes" << std::endl;
    }


    if (false)
    {
        S1.~Setka();
        std::cout << "Setka delete\n";
        return 0;
    }

    S1.Save_cell_parameters("parameters_0209.bin");
    //S1.Save_cell_parameters("parameters_0138.bin");
    //S1.Save_cell_pui_parameters("parameters_0026.bin");

    //S1.Edges_create();
    //S1.Culc_divergence_in_cell();
    //S1.Culc_rotors_in_cell();

    S1.Save_for_interpolate("For_intertpolate_200-.bin", false);
    Interpol SS = Interpol("For_intertpolate_200-.bin");

    //S1.Save_for_interpolate("For_intertpolate_137-.bin", false);
    //Interpol SS = Interpol("For_intertpolate_137-.bin");

    cout << "AAA" << endl;

    if (false)
    {
        // Начальная инициализация
        std::unordered_map<string, double> param;
        std::array<Cell_handle, 6> prev_cell;
        std::array<Cell_handle, 6> next_cell;
        for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();

        SS.Get_param(10.0, 0.0, 0.0, param, prev_cell, next_cell);     // Интерполируем переменные
        for (short int i = 0; i < 6; i++) prev_cell[i] = next_cell[i]; // Обновляем предыдущую ячейку


        cout << "BBB" << endl;
        for (const auto& [key, value] : param) {
            std::cout << key << ":  " << value << '\n';
        }

        cout << "SSSSS" << endl;

        std::unordered_map<string, double> param2;

        SS.Get_HP(10.0, 0.0, 0.0, param2);

        for (const auto& [key, value] : param2) 
        {
            std::cout << key << ":  " << value << '\n';
        }
        return 0.0;

    }

    if (false) // Проверка интерполятора
    {
        // Начальная инициализация
        std::unordered_map<string, double> param;
        std::array<Cell_handle, 6> prev_cell;
        std::array<Cell_handle, 6> next_cell;
        for (short int i = 0; i < 6; i++) prev_cell[i] = Cell_handle();

        std::unordered_map<string, double> parameters;
        // Открываем файл для записи
        std::ofstream outfile("angles.txt");
        if (!outfile.is_open()) {
            return 1;
        }
        const int N = 500;         // Количество шагов
        const double step = const_pi / N;  // Размер шага

        if (false)
        {
            // Основной цикл
            for (double angle = 0.0; angle <= const_pi / 2.0 + 1e-6; angle += step)
            {
                for (double r = 10.0; r < 300.0; r += 0.01)
                {
                    double x = r * cos(angle);
                    double y = r * sin(angle);
                    double z = 0.0;
                    short int zoon = 0;
                    // Записываем в файл
                    SS.Get_param(x, y, z, parameters, zoon);
                    if (zoon == 3)
                    {
                        outfile << x << " " << y << std::endl;
                        break;
                    }
                }
            }

            // Закрываем файл
            outfile.close();
        }

        // Открываем файл для записи
        outfile = std::ofstream("angles2.txt");
        if (!outfile.is_open()) {
            return 1;
        }

        // Основной цикл
        cout << "AA1" << endl;
        for (double angle = 0.0; angle <= const_pi/2; angle += step)
        {
            //cout << "angle = " << angle << endl;
            double x = cos(angle);
            double y = sin(angle);
            double z = 0.0;
            short int zoon = 0;
            // Записываем в файл
            SS.Get_HP(x, y, z, parameters);
            double r = parameters["r"];
            outfile << r * cos(angle) << " " << r * sin(angle) << std::endl;
        }
        cout << "AA3" << endl;

        // Закрываем файл
        outfile.close();

        //std::unordered_map<string, double> param;
        SS.Get_TS(13.94, 0.0, 0.0, param);
        for (const auto& [key, value] : param) {
            std::cout << key << ":  " << value << '\n';
        }
        cout << "A " << endl;


        outfile = std::ofstream("2D.txt");

        for (double x = -200.0; x < 400.0; x = x + 0.4)
        {
            for (double y = -400.0; y < 400.0; y = y + 0.3)
            {
                bool vv = SS.Get_param(x, y, 0.0, parameters, prev_cell, next_cell);
                if (vv == false) continue;
                outfile << x << " " << y << " " << parameters["rho"] << std::endl;
            }
        }

        outfile.close();

        SS.Get_param(0.00001, -350, 0.0, parameters, prev_cell, next_cell);
        for (const auto& [key, value] : parameters) {
            std::cout << key << ":  " << value << '\n';
        }
        cout << "B1 " << endl;

        SS.Get_param(-0.00001, -350, 0.0, parameters, prev_cell, next_cell);
        for (const auto& [key, value] : parameters) {
            std::cout << key << ":  " << value << '\n';
        }
        cout << "B2 " << endl;
        exit(-1);
    }


    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0), "_(1, 0, 0)_", 500.0);

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(-1.0, 0.0, 0.0), "_(-1, 0, 0)_", 500.0);

    S1.Tecplot_print_1D(&SS, Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0), "_(0, 1, 0)_", 500.0);

    S1.Tecplot_print_2D(&SS, 0.0, 0.0, 1.0, -0.00001, "_2d_(0, 0, 1, 0)_");


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


    std::cout << "Hello World!\n";

    S1.~Setka();
    std::cout << "Setka delete\n";


    SS.~Interpol();
    std::cout << "Interpol delete\n";


}

