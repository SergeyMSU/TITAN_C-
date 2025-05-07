// TITAN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "Header.h"
using namespace std;
//class Setka;

int main()
{
    Setka S1 = Setka();

    S1.Read_old_surface("ASurf_Save00591.bin");
    S1.Surf1->Print_TS();
    S1.Surf1->Print_HP();
    S1.Surf1->Print_BS();

    //cout << "TS = " << S1.Surf1->Get_TS(0.0, 0.0) << endl;
    //cout << "HP = " << S1.Surf1->Get_HP(0.0, 0.0, 0) << endl;

    //S1.Move_to_surf(S1.Surf1);
    


    //S1.Tecplot_print_all_yzel_in_3D("SDK1");
    


    //S1.Tecplot_print_all_lush_in_3D("A_Luch");
    //S1.Tecplot_print_all_lush_in_3D("B_Luch");
    //S1.Tecplot_print_all_lush_in_3D("C_Luch");
    //S1.Tecplot_print_all_lush_in_3D("D_Luch");
    //S1.Tecplot_print_all_lush_in_3D("E_Luch");
    //S1.Tecplot_print_all_lush_in_3D("H_Luch");
    //S1.Tecplot_print_all_lush_in_3D("G_Luch");
    //S1.Tecplot_print_all_cell_in_3D();
    //S1.Tecplot_print_krug_yzel_in_3D(2);
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("A_Luch");
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("C_Luch");
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("A2_Luch");
    //S1.Tecplot_print_opor_yzel_in_luchs_3D("C2_Luch");

    //S1.Tecplot_print_all_lush_in_2D();
    std::cout << "Hello World!\n";
}

