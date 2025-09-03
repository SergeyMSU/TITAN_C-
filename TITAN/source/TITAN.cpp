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

    S1.Download_cell_parameters("parameters_0117.bin");   // 107
    // 19 стартовая точка от которой две параллели с пикапами и без
    // 32 с пикапами
    // 62 включи TVD
    // 76 перед тем, как отключить все ТВД и все костыли

    S1.geo->R0 = 0.237455;

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

    //S1.Algoritm(2);
    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");

    

    S1.Find_Yzel_Sosed_for_BS();

    S1.Smooth_angle_HP();
    S1.Smooth_head_HP3();
    S1.Smooth_head_TS3();


    for (int i = 1; i <= 3 * 12; i++) // 6 * 2
    {
        auto start = std::chrono::high_resolution_clock::now();
        cout << "IIIII = " << i << endl;

        S1.Go(false, 400, 1); // 400   1
        S1.Go(true, 200, 1); // 400   1
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

        if (i % 3 == 0)
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

    //S1.Save_cell_parameters("parameters_0052.bin");
    S1.Save_cell_parameters("parameters_0118.bin");
    //S1.Save_cell_pui_parameters("parameters_0026.bin");

    //S1.Edges_create();
    //S1.Culc_divergence_in_cell();
    //S1.Culc_rotors_in_cell();

    S1.Save_for_interpolate("For_intertpolate_71.bin", true);
    Interpol SS = Interpol("For_intertpolate_71.bin");

    if (false) // Проверка интерполятора
    {
        std::unordered_map<string, double> parameters;
        Cell_handle next_cell;
        Cell_handle prev_cell = Cell_handle();
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
                    SS.Get_param(x, y, z, parameters, prev_cell, next_cell, zoon);
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

        std::unordered_map<string, double> param;
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

