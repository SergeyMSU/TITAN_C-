#pragma once
#include "Header.h"

using namespace std;



class Setka
{
public:
	class Geo_param* geo;

	vector<Yzel*> All_Yzel;
	vector<Gran*> All_Gran;
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
	~Setka();
	void New_initial();  // Функция начального построения сетки из файлов 2Д сетки и триангуляции круга
	void New_connect();  // Функция начального построения сетки: создания граней и связывания ячеек
	void Calculating_measure(unsigned short int st_time);
	// Вычисление площадей, центров и нормалей граней, объёмов и центров ячеек

	void Winslow_method(void);  // Реализован winslow метод для триангуляции круга (можно для любой фигуры адаптировать)

	///// БЛОК ГЕОМЕТРИИ
	// Малые вспомогательные функции

	void Renumerate(void); // Перенумерует все элементы сетки (необходимо в случае добавления новых или удаления)
	bool Test_geometr(void);  // Тестируем геометрию сетки (что всё связано и работает корректно)

	void Set_luch_parametr();
	// Добавляет лучам необходимые параметры (например, углы для радиальных лучей) - это ускорит расчёты
	// Нужна для начального построения сетки

	void auto_set_luch_geo_parameter(int for_new);
	// Автоматическая настройка геометрических параметров сетки:
	// настройка для каждого луча (особенно сплайн-лучей) длины начального и конечного вектора
	// работа функции может быть долгой, должна вызываться в программе не часто (не на каждом шаге)

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
	// Печатает лучи в 2Д плоскости, но есть временной параметр, меняя который получаем лучи в разных плоскостях
	// т.е. просматриваем фактически всю сетку (кроме части в головной и хвостовой области)
	void Tecplot_print_plane_lush(int plane); // как предыдущая, но печатает в конкретной плоскости
	void Tecplot_print_plane_surfase(int plane); // как предыдущая, но печатает в конкретной плоскости
	void Tecplot_print_All_surfase_in_2D(); // как предыдущая, но печатает в конкретной плоскости
};

