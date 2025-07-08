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

    S1.Download_cell_parameters("parameters_0006.bin");  // 4
    // c 4 включил bn = 0
    // с 5 начались проблемы с давлением на контакте (становиться меньше 0)
    // с 6 добавил обнуление bn перед контактом

    cout << "C2 " << endl;

    S1.auto_set_luch_geo_parameter(0);

    S1.Init_TVD();
    cout << "D2 " << endl;

    S1.Init_physics();

    cout << "E " << endl;

    //S1.Smooth_head_TS2();
    ///S1.Smooth_head_HP2();

    S1.Tecplot_print_cell_plane_parameters();
    S1.Tecplot_print_all_lush_in_2D();
    S1.Tecplot_print_all_cell_in_3D();

    //S1.Algoritm(2);
    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");
    

    for (int i = 1; i <= 6 * 11; i++) // 6 * 2
    {
        auto start = std::chrono::high_resolution_clock::now();
        cout << "IIIII = " << i << endl;
        S1.Go(false, 400, 1); // 400

        S1.Tecplot_print_cell_plane_parameters();
        S1.Tecplot_print_all_lush_in_2D();
        S1.Tecplot_print_all_gran_in_surface("TS");
        S1.Tecplot_print_all_gran_in_surface("HP");
        S1.Tecplot_print_all_gran_in_surface("BS");
        //S1.Go(true, 100, 1);
        //S1.Tecplot_print_cell_plane_parameters();

        //S1.Init_physics();

        if (i % 6 == 0)
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

    S1.Save_cell_parameters("parameters_0007.bin");

    S1.Save_for_interpolate("For_intertpolate_1.bin");
    Interpol SS = Interpol("For_intertpolate_1.bin");

    

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

