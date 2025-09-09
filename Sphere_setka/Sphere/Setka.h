#pragma once
using namespace std;
#include "Header.h"

class point;
class Cell;


class Setka
{
	double phi = pi * 30.0 / 180.0;
	int Nphi = 60;     // Сколько точек по кругу (вращений сетки вокруг оси Х)
	int Ntheta = 40;        // Сколько точек по другому углу в добавочном слое
public:
	vector<point*> all_points;
	vector<Cell*> all_Cells;
	point*** PP;   // Часть точек, которая в сферических координатах (она дублирует точки в общем массиве, но удобна для 
	// быстрого и понятного доступа)

	void Intitial_read(void);
	void regularize(void); // Регуляризация точек в кругу
	void regularize2(void); // Регуляризация точек в кругу
	void Intitial_build(void);  // Проецируем точки на сфере и создаём точки с другой стороны
	void Intitial_build2(void); // Создаём остальные точки на сфере
	void Intitial_build3(void); // Создаём остальные точки на сфере

	void Find_cell_soseds(void); // Находим соседей ячейкам
	bool Cell_is_soseds(Cell* A, Cell* B);

	point* Find_point_id(const int& num);
	point* Get_point(const double& x, const double& y, const double& z, const double& r);

	void Print_Setka_TecPlot(void);
	void Print_Setka_TecPlot_surface(void);
	void Print_cell_soseds(void);

	void All_numerate(); // Нумерация узлов и ячеек
	void Save_for_3D(string name);
};

