#pragma once
#include "Header.h"

typedef KKexact::Point_3 Point;
typedef KKexact::Vector_3 Vector;
typedef KKexact::Tetrahedron_3 Tetrahedron;
typedef KKexact::Point_2 Point2;

class Interpol
{
public:
	std::vector<std::pair<Point, size_t>> points;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point2, size_t>> point_TS;  // Точки и их номера для триангуляции
	std::vector <Int_point*> Cells;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_TS;     // Точки в которых хранятся параметры
	Delaunay* Delone;
	Delaunay2* Delone_TS;

	vector<string> param_names;  // Названия всех хранящихся переменных

	unordered_map<string, double> stepen;   // Стемени радиуса в интерполяции

	double L6;  // До какого расстояния слева выделяется контакт

	Interpol(string name);
	~Interpol();


	bool Get_param(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters, const Cell_handle& prev_cell, Cell_handle& next_cell);

	bool Get_param(const double& x, const double& y, const double& z, 
		std::unordered_map<string, double>& parameters, 
		const Cell_handle& prev_cell, Cell_handle& next_cell, short int& this_zone);

	bool Get_TS(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
};

