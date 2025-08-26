#pragma once
#include "Header.h"

typedef KKexact::Point_3 Point;
typedef KKexact::Vector_3 Vector;
typedef KKexact::Tetrahedron_3 Tetrahedron;
typedef KKexact::Point_2 Point2;



class Interpol
{
public:
	std::vector<std::pair<Point, size_t>> points_1;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point, size_t>> points_2;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point, size_t>> points_3;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point, size_t>> points_4;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point, size_t>> points_5;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point, size_t>> points_6;  // ����� � �� ������ ��� ������������

	std::vector<std::pair<Point2, size_t>> point_TS;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point2, size_t>> point_BS;  // ����� � �� ������ ��� ������������
	std::vector<std::pair<Point2, size_t>> point_HP_1;  // ����� � �� ������ ��� ������������

	std::vector <Int_point*> Cells_1;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_2;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_3;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_4;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_5;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_6;     // ����� � ������� �������� ���������

	std::vector <Int_point*> Cells_TS;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_BS;     // ����� � ������� �������� ���������
	std::vector <Int_point*> Cells_HP_1;     // ����� � ������� �������� ���������
	boost::multi_array<Int_point*, 2> Cells_HP_2;


	Delaunay* Delone_1;
	Delaunay* Delone_2;
	Delaunay* Delone_3;
	Delaunay* Delone_4;
	Delaunay* Delone_5;
	Delaunay* Delone_6;

	Delaunay2* Delone_TS;
	Delaunay2* Delone_BS;
	Delaunay2* Delone_HP_1;

	vector<string> param_names;  // �������� ���� ���������� ����������

	unordered_map<string, double> stepen;   // ������� ������� � ������������

	double L6;  // �� ������ ���������� ����� ���������� �������

	Interpol(string name);
	~Interpol();


	bool Get_param(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters, const Cell_handle& prev_cell, Cell_handle& next_cell);

	bool Get_param(const double& x, const double& y, const double& z, 
		std::unordered_map<string, double>& parameters, 
		const Cell_handle& prev_cell, Cell_handle& next_cell, short int& this_zone);

	bool Get_TS(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
	bool Get_HP(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
	bool Get_BS(const double& x, const double& y, const double& z,
		std::unordered_map<string, double>& parameters);
};

