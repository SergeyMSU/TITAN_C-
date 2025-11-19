#include "Setka.h"


bool Setka::Get_pui_SS(vector<double>& pui_Sm, vector<double>& pui_Sp1, vector<double>& pui_Sp2, 
	short int ii, const double& x, const double& y, const double& z,
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell(4);
	vector<double> koeff_cell(4);

	//cout << "Get_pui_Sm  1" << endl;
	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);
	//cout << "Get_pui_Sm  2" << endl;

	if (b == false) return false;

	prev_cell = next_cell;

	Cell* CC;

	for (int jj = 0; jj < this->phys_param->pui_nW; jj++)
	{
		double S = 0.0;
		double SS = 0.0;
		double SSS = 0.0;

		for (int i = 0; i < num_cell.size(); i++)
		{
			int j = num_cell[i];
			double k = koeff_cell[i];

			CC = S_MK.All_Cell[j];
			S += CC->pui_Sm[jj] * k;
			SS += CC->pui_Sp(0, jj) * k;
			if (ii == 2)
			{
				int rows = CC->pui_Sp.rows();
				if (rows >= 2)
				{
					SSS += CC->pui_Sp(1, jj) * k;
				}
				else
				{
					SSS += 0.0;
				}
			}
		}

		pui_Sm[jj] = S;
		pui_Sp1[jj] = SS;
		if (ii == 2) pui_Sp2[jj] = SSS;
	}

	return true;
}

bool Setka::Get_pui_Sm(double& pui_Sm, int n, const double& x, const double& y, const double& z,
	Setka& S_MK, Interpol& SI_MK, Cell_handle& prev_cell, Cell_handle& next_cell)
{
	vector<int> num_cell(4);
	vector<double> koeff_cell(4);

	//cout << "Get_pui_Sm  1" << endl;
	bool b = SI_MK.Get_real_cells(x, y, z, num_cell, koeff_cell, prev_cell, next_cell);
	//cout << "Get_pui_Sm  2" << endl;

	if (b == false) return false;

	prev_cell = next_cell;

	double S = 0.0;
	//cout << "Get_pui_Sm  3" << endl;
	for (int i = 0; i < num_cell.size(); i++)
	{
		int j = num_cell[i];
		double k = koeff_cell[i];

		//cout << "Get_pui_Sm  4" << endl;

		//cout << "j = " << j << "   " << S_MK.All_Cell.size() << endl;
		//cout << S_MK.All_Cell[j]->pui_Sm.size() << "   " << n << endl;
		S += S_MK.All_Cell[j]->pui_Sm[n] * k;


		//cout << "Get_pui_Sm  5" << endl;
	}
	//cout << "Get_pui_Sm  6" << endl;

	pui_Sm = S;
	return true;
}


bool Setka::Get_pui_Sp(double& pui_Sp, short int ii, int n, const double& x, const double& y, const double& z,
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