// TITAN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "Header.h"
using namespace std;
//class Setka;

int main()
{
    cout << "Start Programm" << endl;
    Setka S1 = Setka();

    //S1.phys_param->raspad_testing();

    //return 0;

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
    S1.Init_physics();
    cout << "D " << endl;
    S1.Init_TVD();
    cout << "D2 " << endl;

    S1.Download_cell_parameters("parameters_0003.bin");

    cout << "E " << endl;

    //S1.Go(true, 300, 1);
    //S1.Go(false, 300, 1);
    //S1.Tecplot_print_cell_plane_parameters();
    for (int i = 0; i < 10; i++)
    {
        cout << "IIIII = " << i << endl;
        S1.Go(false, 10000, 1);
        S1.Go(true, 300, 1);
        S1.Tecplot_print_cell_plane_parameters();
        S1.Save_cell_parameters("parameters_0000.bin");
    }




    S1.Save_cell_parameters("parameters_0004.bin");

    cout << "F " << endl;

    S1.Tecplot_print_cell_plane_parameters();


    cout << "YSPEX" << endl;


    S1.Tecplot_print_all_yzel_in_3D("SDK1");
    S1.Tecplot_print_gran_with_condition();
    


    //S1.Tecplot_print_all_lush_in_3D("A_Luch");
    //S1.Tecplot_print_all_lush_in_3D("B_Luch");
    //S1.Tecplot_print_all_lush_in_3D("C_Luch");
    //S1.Tecplot_print_all_lush_in_3D("D_Luch");
    //S1.Tecplot_print_all_lush_in_3D("E_Luch");
    //S1.Tecplot_print_all_lush_in_3D("H_Luch");
    //S1.Tecplot_print_all_lush_in_3D("G_Luch");
    S1.Tecplot_print_all_cell_in_3D();
    S1.Tecplot_print_krug_yzel_in_3D(1);
    S1.Tecplot_print_krug_yzel_in_3D(2);
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("A_Luch");
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("C_Luch");
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("A2_Luch");
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("C2_Luch");

    S1.Tecplot_print_all_lush_in_2D();
    S1.Tecplot_print_All_surfase_in_2D();
    S1.Tecplot_print_plane_lush(0);
    S1.Tecplot_print_plane_surfase(0);
    S1.Tecplot_print_all_gran_in_cell();
    S1.Tecplot_print_all_gran_in_surface("TS");
    S1.Tecplot_print_all_gran_in_surface("HP");
    S1.Tecplot_print_all_gran_in_surface("BS");
    S1.Tecplot_print_all_yzel_with_condition();


    std::cout << "Hello World!\n";
}

