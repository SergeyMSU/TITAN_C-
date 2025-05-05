#pragma once
#include "Header.h"

using namespace std;

class Setka
{
public:
	class Geo_param* geo;

	vector<Yzel*> All_Yzel;
	vector <vector<Yzel*>> Yzel_2D;
	vector <vector<Yzel*>> Krug_Yzel;         // Узлы в кругу (в головной области)  
	// первая координата - на какой сфере точке (0 - для ближник точек), вторая координата - номер точки на луче
	// Krug_Yzel и Yzel_2D пересекаются (крайтими узлами). Однако A_Luch и A2_Luch не пересекаются
	vector <vector<Yzel*>> Krug_Yzel_2;       // Узлы в кругу (в хвосте)

	vector<Cell*> All_Cell;           // Все ячейки сетки
	vector <vector<Cell*>> Cell_2D;   // Ячейки в части сетки, которая получена вращением 2Д сетки
	vector <vector<Cell*>> Cell_layer_head;   // Ячейки в головной части сетки (вблизи оси)
	vector <vector<Cell*>> Cell_layer_tail;   // Ячейки в хвостовой части сетки (вблизи оси)

	vector<Luch*> All_Luch;
	// Лучи для конкретной реализации сетки (не универсальный блок программы)
	vector <vector<Luch*>> A_Luch;   //  [Plane][Number]  сначала номер плоскости, в которой лежит луч, затем номер самого луча
	vector <vector<Luch*>> B_Luch;
	vector <vector<Luch*>> C_Luch;
	vector <vector<Luch*>> D_Luch;
	vector <vector<Luch*>> E_Luch;
	vector <vector<Luch*>> H_Luch;
	vector <vector<Luch*>> G_Luch;
	map <string, vector<vector<Luch*>>*> All_name_luch; 
	vector<string> name_luch;

	vector<Luch*> A2_Luch;
	vector<Luch*> C2_Luch;


	Surfaces* Surf1;
	

	Setka();
	void New_initial();  // Функция начального построения сетки из файлов 2Д сетки и триангуляции круга

	void Winslow_method(void);  // Реализован winslow метод для триангуляции круга (можно для любой фигуры адаптировать)

	///// БЛОК ГЕОМЕТРИИ
	// Малые вспомогательные функции

	void Renumerate(void); // Перенумерует все элементы сетки (необходимо в случае добавления новых или удаления)

	void Set_luch_parametr();
	// Добавляет лучам необходимые параметры (например, углы для радиальных лучей) - это ускорит расчёты

	void Read_old_surface(string name);
	// Считывает файл поверхностей (в формате старой программы на фортране)

	void Move_to_surf(Surfaces* Surf);

	// ВИЗУАЛИЗАЦИЯ
	// Ниже всё, что касается визуализации сетки

	// Для Tecplot
	void Tecplot_print_all_yzel_in_3D(string name);
	// Печатает все узлы (но по слоям, чтобы можно было удобно их просмотреть)

	void Tecplot_print_all_lush_in_3D(string name);
	void Tecplot_print_all_cell_in_3D();
	void Tecplot_print_krug_yzel_in_3D(int num);
	void Tecplot_print_opor_yzel_in_luchs_3D(string name);
	void Tecplot_print_all_lush_in_2D();
};

