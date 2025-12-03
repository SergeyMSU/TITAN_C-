// TITAN.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#include "Header.h"
using namespace std;
//class Setka;

int main()
{
    cout << "Start Programm" << endl;

    // Создаём основную сетку из файлов вспомогательных сеток
    Setka S1 = Setka("SDK1_2D_Setka.bin", "SDK1_krug_setka.bin", 60);

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

    // Считываем физические параметры и геометрическое положение узлов из файла (предыдущего расчёта)
    //S1.Download_cell_parameters("parameters_0060.bin");   

    S1.Download_cell_parameters("parameters_0065.bin");


    //S1.Download_cell_parameters("parameters_promeg_1124.bin");

    // Ещё один блок обязательной настройки
    if (true)
    {
        // Точно задаём положение внутренней границы сетки
        S1.geo->R0 = S1.phys_param->R_0;

        // Автоматически подстраиваем геометрические параметры сетки под новые положения узлов
        S1.auto_set_luch_geo_parameter(0);

        // Инициализируем TVD (находим соседей и т.д.)
        S1.Init_TVD();
    }


    //  Ручное изменение BS
    if (false)
    {
        cout << "Hand" << endl;
        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        S1.Culc_Velocity_surface(0, 0.0, 3);
        for (int i_step = 0; i_step < S1.All_Luch.size(); i_step++)
        {
            auto lu = S1.All_Luch[i_step];
            lu->dvigenie(1);
        }
        for (auto& i : S1.All_Yzel)
        {
            i->coord[0][0] = i->coord[1][0];
            i->coord[0][1] = i->coord[1][1];
            i->coord[0][2] = i->coord[1][2];
        }
        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        S1.auto_set_luch_geo_parameter(0);
        for (auto& i : S1.All_Yzel)
        {
            i->coord[1][0] = i->coord[0][0];
            i->coord[1][1] = i->coord[0][1];
            i->coord[1][2] = i->coord[0][2];
        }
        S1.Calculating_measure(0);
        S1.Calculating_measure(1);
        S1.Init_TVD();
    }

    cout << "A-" << endl;
    // Задаём начальные и граничные условия
    S1.Init_physics();

    // Блок начальной визуализации сетки для проверки корректного построения
    if (true)
    {
        S1.Tecplot_print_cell_plane_parameters();
        S1.Tecplot_print_all_lush_in_2D();
        S1.Tecplot_print_2D_setka(0.0, 0.0, 1.0, -0.00001, "init_setka_2d_(0, 0, 1, 0)_");
        S1.Tecplot_print_2D_setka(0.0, 1.0, 0.0, -0.00001, "init_setka_2d_(0, 1, 0, 0)_");
        S1.Tecplot_print_2D_setka(0.0, 1.0, 1.0, -0.00001, "init_setka_2d_(0, 1, 1, 0)_");
        S1.Tecplot_print_all_gran_in_surface("TS");
        S1.Tecplot_print_all_gran_in_surface("HP");
        S1.Tecplot_print_all_gran_in_surface("BS");
    }

    // Выбор основного алгоритма расчёта (в данной функции представлены все варианты расчёта: атомы, мгд и т.д.), см. саму функцию
    S1.Algoritm(10, &S1);
    //S1.Algoritm(8, &S1);

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



    //S1.Save_cell_parameters("parameters_0065.bin");
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

    //S1.Save_for_interpolate("For_intertpolate_0065.bin", true);
    //Interpol SS = Interpol("For_intertpolate_0064.bin");

    S1.Save_for_interpolate("For_intertpolate_0059-.bin", false);
    Interpol SS = Interpol("For_intertpolate_0059-.bin");

    cout << "AAA" << endl;

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
}

