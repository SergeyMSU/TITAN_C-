#pragma once
#include "Header.h"

typedef KKexact::Point_3 Point;
typedef KKexact::Vector_3 Vector;
typedef KKexact::Tetrahedron_3 Tetrahedron;

class Interpol
{
public:
	std::vector<std::pair<Point, size_t>> points;  // ����� � �� ������ ��� ������������
	std::vector <Int_point*> Cells;     // ����� � ������� �������� ���������
	Delaunay* Delone;
	vector<string> param_names;  // �������� ���� ���������� ����������


	Interpol(string name);
	~Interpol();

	bool Get_param(const double& x, const double& y, const double& z, 
		std::unordered_map<string, double>& parameters, 
		const Cell_handle& prev_cell, Cell_handle& next_cell);
};

