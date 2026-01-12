#pragma once
#include "Header.h"

//typedef KKexact::Point_3 Point;
//typedef KKexact::Vector_3 Vector;
//typedef KKexact::Tetrahedron_3 Tetrahedron;
//typedef KKexact::Point_2 Point2;



class Interpol
{
public:
	std::vector<std::pair<Point, size_t>> points_1;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point, size_t>> points_2;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point, size_t>> points_3;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point, size_t>> points_4;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point, size_t>> points_5;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point, size_t>> points_6;  // Точки и их номера для триангуляции

	std::vector<std::pair<Point2, size_t>> point_TS;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point2, size_t>> point_BS;  // Точки и их номера для триангуляции
	std::vector<std::pair<Point2, size_t>> point_HP_1;  // Точки и их номера для триангуляции

	std::vector <Int_point*> Cells_1;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_2;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_3;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_4;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_5;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_6;     // Точки в которых хранятся параметры

	std::vector <Int_point*> Cells_TS;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_BS;     // Точки в которых хранятся параметры
	std::vector <Int_point*> Cells_HP_1;     // Точки в которых хранятся параметры
	boost::multi_array<Int_point*, 2> Cells_HP_2;

	mutex mut_Delone_1;
	Delaunay* Delone_1 = nullptr;
	Delaunay* Delone_2 = nullptr;
	Delaunay* Delone_3 = nullptr;
	Delaunay* Delone_4 = nullptr;
	Delaunay* Delone_5 = nullptr;
	Delaunay* Delone_6 = nullptr;

	Delaunay2* Delone_TS = nullptr;
	Delaunay2* Delone_BS = nullptr;
	Delaunay2* Delone_HP_1 = nullptr;

	vector<string> param_names;  // Названия всех хранящихся переменных

	unordered_map<string, double> stepen;   // Степени радиуса в интерполяции

	double L6;  // До какого расстояния слева выделяется контакт
	bool razriv;       // Есть ли особая интерполяция на разрывах?

	Interpol(string name);
	~Interpol();


	bool Get_param(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
	bool Get_param(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters, const std::array<Cell_handle, 6>& prev_cell,
		std::array<Cell_handle, 6>& next_cell,
		short int& this_zone);
	bool Get_param(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters, short int& this_zone);
	bool Get_param(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters, const std::array<Cell_handle, 6>& prev_cell,
		std::array<Cell_handle, 6>& next_cell);

	bool Get_TS(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
	bool Get_HP(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
	bool Get_BS(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);

	bool Get_real_cells(const double& x, const double& y, const double& z,
		vector<int>& num_cell, vector<double>& koeff_cell, const Cell_handle& prev_cell, Cell_handle& next_cell);
	// работает только для razriv = false 
	// Для данной точки возвращает набор из 4-х номеров ячеек в исходной сетке
	// и набор из 4-х коэффициентов для ручной интерполяции нужных переменных 
};

