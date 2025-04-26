#pragma once
#include <vector>
#include "Yzel.h"
#include "Header.h"

using namespace std;


class Setka
{
public:
	struct Geo_param geo;           // �������������� ��������� �����



	vector<Yzel*> All_Yzel;
	vector<Luch*> All_Luch;
	vector<Luch*> A_Luch;
	vector<Luch*> B_Luch;
	vector<Luch*> C_Luch;
	vector<Luch*> D_Luch;
	vector<Luch*> E_Luch;
	vector<Luch*> H_Luch;
	vector<Luch*> G_Luch;
	vector < vector<Luch*>*> All_name_luch;

	vector<Cell*> All_Cell;

	vector<vector<Cell*>> A_Cells;
	vector<vector<Cell*>> B_Cells;
	vector<vector<Cell*>> C_Cells;
	vector<vector<Cell*>> D_Cells;
	vector<vector<Cell*>> E_Cells;
	vector<vector<Cell*>> H_Cells;
	vector<vector<Cell*>> G_Cells;

	Setka();
	void Set_geo();
	void Construct_initial();
	void All_numerate(); // ��������� ����� � �����

	bool Cell_is_soseds(Cell* A, Cell* B);  // ���������, �������� �� ������ ������ ��������

	void Print_yzel();        // ������ ���� �����
	void Print_yzels_opor(string name);    // ������ ���� ������� ����� ��� ���� name
	void Print_cell();    // ������ ���� ����� 
	void Print_luch();    // ������ ���� �����
	void Print_cell_soseds();    // ������ ���� �����
	void Print_cell_center();    // ������ ������� �����
	void Print_cell_center_test();    // ��� ������������ ���������, ������ ������� ����� � �����-�� ��������

	void Save_for_3D(string name); // ���������� ����� ��� 3� ���������
	void Testing(void); // ����� ������������ ����� ... ����� ����� ��� �����

};

