#include "Setka.h"

bool Setka::Get_pui_Sm(double& pui_Sm, int n, double& x, double& y, double& z, 
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell;
	vector<double> koeff_cell;

	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);

	if (b == false) return false;

	prev_cell = next_cell;

	double S = 0.0;
	for (int i = 0; i < num_cell.size(); i++)
	{
		int j = num_cell[i];
		double k = koeff_cell[i];

		S += S_MK.All_Cell[j]->pui_Sm[n] * k;
	}

	pui_Sm = S;
	return true;
}


bool Setka::Get_pui_Sp(double& pui_Sp, short int ii, int n, double& x, double& y, double& z,
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell;
	vector<double> koeff_cell;

	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);

	if (b == false) return false;

	prev_cell = next_cell;

	double S = 0.0;
	for (int i = 0; i < num_cell.size(); i++)
	{
		int j = num_cell[i];
		double k = koeff_cell[i];

		S += S_MK.All_Cell[j]->pui_Sp(ii, n) * k;
	}

	pui_Sp = S;
	return true;
}