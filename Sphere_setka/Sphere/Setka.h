#pragma once
using namespace std;
#include "Header.h"

class point;
class Cell;


class Setka
{
	double phi = pi * 30.0 / 180.0;
	int Nphi = 60;     // ������� ����� �� ����� (�������� ����� ������ ��� �)
	int Ntheta = 40;        // ������� ����� �� ������� ���� � ���������� ����
public:
	vector<point*> all_points;
	vector<Cell*> all_Cells;
	point*** PP;   // ����� �����, ������� � ����������� ����������� (��� ��������� ����� � ����� �������, �� ������ ��� 
	// �������� � ��������� �������)

	void Intitial_read(void);
	void regularize(void); // ������������� ����� � �����
	void regularize2(void); // ������������� ����� � �����
	void Intitial_build(void);  // ���������� ����� �� ����� � ������ ����� � ������ �������
	void Intitial_build2(void); // ������ ��������� ����� �� �����
	void Intitial_build3(void); // ������ ��������� ����� �� �����

	void Find_cell_soseds(void); // ������� ������� �������
	bool Cell_is_soseds(Cell* A, Cell* B);

	point* Find_point_id(const int& num);
	point* Get_point(const double& x, const double& y, const double& z, const double& r);

	void Print_Setka_TecPlot(void);
	void Print_Setka_TecPlot_surface(void);
	void Print_cell_soseds(void);

	void All_numerate(); // ��������� ����� � �����
	void Save_for_3D(string name);
};

