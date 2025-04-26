// 2D_setka_construct.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Header.h"

int main()
{
    Setka S = Setka();

    S.Set_geo();
    S.Construct_initial();
    S.Print_yzel();
    S.Print_cell();
    S.Print_cell_center();
    S.Print_luch();
    S.Print_yzels_opor("A");
    S.Print_yzels_opor("B");
    S.Print_yzels_opor("C");
    S.Print_yzels_opor("D");
    S.Print_yzels_opor("E");
    S.Print_yzels_opor("H");
    S.Print_yzels_opor("G");
    S.Print_cell_soseds();
    S.Print_cell_center_test();
    S.Save_for_3D("SDK1");

    cout << "All cell = " << S.All_Cell.size() << endl;


}
