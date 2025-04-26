#pragma once
#include <vector>
#include "Yzel.h"
#include "Header.h"

using namespace std;


class Setka
{
public:
	struct Geo_param geo;           // √еометрические параметры сетки



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
	void All_numerate(); // Ќумераци€ узлов и €чеек

	bool Cell_is_soseds(Cell* A, Cell* B);  // ѕроверить, €вл€ютс€ ли данные €чейки сосед€ми

	void Print_yzel();        // ѕечать всех узлов
	void Print_yzels_opor(string name);    // ѕечать всех опорных узлов дл€ луча name
	void Print_cell();    // ѕечать всех €чеек 
	void Print_luch();    // ѕечать всех лучей
	void Print_cell_soseds();    // ѕечать всех лучей
	void Print_cell_center();    // ѕечать центров €чеек
	void Print_cell_center_test();    // ƒл€ тестировани€ геометрии, ѕечать центров €чеек с каким-то условием

	void Save_for_3D(string name); // —охранение сетки дл€ 3ƒ программы
	void Testing(void); // ќбщее тестирование сетки ... здесь будут все тесты

};

