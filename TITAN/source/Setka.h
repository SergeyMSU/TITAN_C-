#pragma once
#include "Header.h"

using namespace std;

class Setka
{
public:
	class Geo_param* geo;

	vector<Yzel*> All_Yzel;
	vector <vector<Yzel*>> Yzel_2D;
	vector <vector<Yzel*>> Krug_Yzel;         // ���� � ����� (� �������� �������)  
	// ������ ���������� - �� ����� ����� ����� (0 - ��� ������� �����), ������ ���������� - ����� ����� �� ����
	// Krug_Yzel � Yzel_2D ������������ (�������� ������). ������ A_Luch � A2_Luch �� ������������
	vector <vector<Yzel*>> Krug_Yzel_2;       // ���� � ����� (� ������)

	vector<Cell*> All_Cell;           // ��� ������ �����
	vector <vector<Cell*>> Cell_2D;   // ������ � ����� �����, ������� �������� ��������� 2� �����
	vector <vector<Cell*>> Cell_layer_head;   // ������ � �������� ����� ����� (������ ���)
	vector <vector<Cell*>> Cell_layer_tail;   // ������ � ��������� ����� ����� (������ ���)

	vector<Luch*> All_Luch;
	// ���� ��� ���������� ���������� ����� (�� ������������� ���� ���������)
	vector <vector<Luch*>> A_Luch;   //  [Plane][Number]  ������� ����� ���������, � ������� ����� ���, ����� ����� ������ ����
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
	void New_initial();  // ������� ���������� ���������� ����� �� ������ 2� ����� � ������������ �����

	void Winslow_method(void);  // ���������� winslow ����� ��� ������������ ����� (����� ��� ����� ������ ������������)

	///// ���� ���������
	// ����� ��������������� �������

	void Renumerate(void); // ������������ ��� �������� ����� (���������� � ������ ���������� ����� ��� ��������)

	void Set_luch_parametr();
	// ��������� ����� ����������� ��������� (��������, ���� ��� ���������� �����) - ��� ������� �������

	void Read_old_surface(string name);
	// ��������� ���� ������������ (� ������� ������ ��������� �� ��������)

	void Move_to_surf(Surfaces* Surf);

	// ������������
	// ���� ��, ��� �������� ������������ �����

	// ��� Tecplot
	void Tecplot_print_all_yzel_in_3D(string name);
	// �������� ��� ���� (�� �� �����, ����� ����� ���� ������ �� �����������)

	void Tecplot_print_all_lush_in_3D(string name);
	void Tecplot_print_all_cell_in_3D();
	void Tecplot_print_krug_yzel_in_3D(int num);
	void Tecplot_print_opor_yzel_in_luchs_3D(string name);
	void Tecplot_print_all_lush_in_2D();
};

