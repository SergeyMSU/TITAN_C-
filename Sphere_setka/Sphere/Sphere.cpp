// Sphere.cpp : Работа с построение сетки на сфере для последующего использования в программе построения основной сетки


#include "Header.h"


using namespace std;

int main()
{
    Setka S = Setka();
    S.Intitial_read();

    //A->cells[4]->points[2] = C;
    //C->cells.push_back(A->cells[4]);


    S.regularize();

    // Этот блок для создания 3Д Геометрии
    //S.Intitial_build();
    //S.Intitial_build2();
    //S.Intitial_build3();

    S.Find_cell_soseds();

    S.Print_Setka_TecPlot();
    S.Print_Setka_TecPlot_surface();
    S.Print_cell_soseds();
    S.Save_for_3D("SDK1");
}
